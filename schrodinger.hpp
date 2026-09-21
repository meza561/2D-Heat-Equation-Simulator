#ifndef SCHRODINGER_HPP
#define SCHRODINGER_HPP

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

#include "tridiagonal.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * 2D time-dependent Schrodinger equation, in natural units (hbar = m = 1).
 *
 *     i dpsi/dt = -1/2 lap(psi) + V(x,y) psi
 *
 * The structural difference from every other solver here is that the state is
 * COMPLEX. Heat, reaction-diffusion and wave all evolve a real field; this one
 * evolves psi = u + iv, and the observable is |psi|^2. The field is stored as
 * two parallel real grids rather than a grid of std::complex so that the
 * export path, the database schema and the frontend keep taking exactly the
 * real matrices they already take. Complex arithmetic appears only inside the
 * sweeps, as scratch, where it is genuinely easier to read than hand-expanded
 * real and imaginary parts.
 *
 * WHY CRANK-NICOLSON, AND NOT AN EXPLICIT SCHEME
 *
 * The evolution operator exp(-iHt) is unitary: total probability is conserved
 * exactly, forever. That is not a nicety, it is the physics - |psi|^2 is a
 * probability density, and a scheme that leaks it is computing something that
 * is not quantum mechanics.
 *
 * Forward Euler applied to this equation has amplification factor |1 - iEdt|
 * = sqrt(1 + E^2 dt^2) > 1 for every mode, so it is unconditionally UNSTABLE:
 * the norm grows without bound, fastest for the highest-energy modes. That is
 * the opposite of the heat equation, where explicit stepping is merely
 * conditionally stable and its errors are dissipative - they damp toward the
 * steady state the solution was heading for anyway, so they hide. Here any
 * norm error is immediately visible as probability appearing from nowhere.
 *
 * Crank-Nicolson instead uses the Cayley transform
 *
 *     psi^(n+1) = (I + i dt/2 H)^-1 (I - i dt/2 H) psi^n
 *
 * which for Hermitian H is exactly unitary in exact arithmetic - not
 * approximately, not to O(dt^2). The norm is then conserved to rounding error,
 * which is what the norm-conservation test in tests.cpp pins down.
 *
 * ADI SPLITTING, AND WHY IT COSTS NOTHING HERE
 *
 * A direct 2D Crank-Nicolson solve would need a full sparse matrix. Instead H
 * splits as Hx + Hy and each direction is advanced by its own Cayley
 * transform, one tridiagonal system per grid line, exactly as the heat solver
 * does for its ADI sweeps - the Thomas solver itself is shared, see
 * tridiagonal.hpp. Each half is separately unitary, so their product is
 * unitary too, and the norm survives the splitting intact.
 *
 * Splitting normally costs accuracy: exp(-iHx dt) exp(-iHy dt) differs from
 * exp(-iH dt) by the commutator [Hx, Hy]. Here it costs nothing, because every
 * potential offered is SEPARABLE - V(x,y) = Vx(x) + Vy(y) - so Hx acts only on
 * x, Hy only on y, and the two commute exactly. The scheme is then second
 * order in both dt and dx with no splitting error at all. Hence the potential
 * is stored as two 1D arrays rather than a 2D grid: that representation makes
 * the separability an invariant of the type rather than an assumption.
 * A non-separable potential (a barrier with a slit, say) would need the V/2
 * split across both sweeps, reintroducing an O(dt^3) commutator term.
 *
 * BOUNDARIES
 *
 * psi = 0 on all four edges, for every preset - the domain is a box with hard
 * walls. This is what makes the discrete Hamiltonian Hermitian, and therefore
 * the Cayley transform unitary. For the square well those walls ARE the
 * physics. For the free particle and the harmonic oscillator the box is meant
 * to be large enough that the packet never reaches them; a packet that does
 * reach a wall will reflect off it, which is a real solution of this problem
 * but not of the infinite-domain one the closed forms describe.
 */
class SchrodingerSimulation {
public:
    /**
     * Free         V = 0. A packet spreads and, given long enough, reflects.
     * SquareWell   V = 0 inside; the domain boundary is the infinite wall, so
     *              this shares its potential with Free and differs in which
     *              initial state is natural - a standing eigenstate, not a
     *              travelling packet.
     * Harmonic     V = 1/2 omega^2 r^2 about the box centre.
     */
    enum class Potential { Free, SquareWell, Harmonic };

    using Grid = std::vector<std::vector<double>>;
    using Complex = std::complex<double>;

    SchrodingerSimulation(int rows, int cols, double dx, double dt)
        : rows(rows), cols(cols), dx(dx), dt(dt) {
        u.assign(rows, std::vector<double>(cols, 0.0));
        v = u;
        prob = u;
        vx.assign(cols, 0.0);
        vy.assign(rows, 0.0);
    }

    /**
     * Select the potential. omega is the oscillator frequency and is ignored
     * by the other two presets. Stored as the separable pair Vx(x), Vy(y);
     * see the note on commuting operators above.
     */
    void setPotential(Potential p, double omegaIn = 0.05) {
        pot = p;
        omega = (omegaIn > 0.0) ? omegaIn : 0.05;
        std::fill(vx.begin(), vx.end(), 0.0);
        std::fill(vy.begin(), vy.end(), 0.0);
        if (p == Potential::Harmonic) {
            const double half = 0.5 * omega * omega;
            for (int j = 0; j < cols; ++j) vx[j] = half * xAt(j) * xAt(j);
            for (int i = 0; i < rows; ++i) vy[i] = half * yAt(i) * yAt(i);
        }
        // Free and SquareWell both leave V at zero: the well's walls are the
        // Dirichlet boundary every preset already has.
    }

    Potential potential() const { return pot; }

    /** Physical coordinates. The origin sits at the centre of the box, so a
     *  packet seeds at (0, 0) and the harmonic well is 1/2 omega^2 r^2 about
     *  it, rather than every caller having to offset by half the grid. */
    double xAt(int j) const { return (j - 0.5 * (cols - 1)) * dx; }
    double yAt(int i) const { return (i - 0.5 * (rows - 1)) * dx; }

    int getRows() const { return rows; }
    int getCols() const { return cols; }

    /** Width of the oscillator ground state, sigma = 1/sqrt(2 omega). Only
     *  meaningful for the Harmonic preset - the other two have no bound
     *  state of their own. */
    double groundStateSigma() const { return 1.0 / std::sqrt(2.0 * omega); }

    /**
     * A Gaussian wave packet, normalised to unit total probability.
     *
     *     psi = exp(-r^2 / (4 sigma0^2)) exp(i (kx x + ky y))
     *
     * The 4 in the exponent is the convention that makes sigma0 the standard
     * deviation of the DENSITY |psi|^2 rather than of the amplitude, which is
     * what the free-particle spreading law sigma(t)^2 = sigma0^2 + (t/2sigma0)^2
     * is written in. (kx, ky) is the mean momentum, so the packet drifts at
     * that velocity; it is zero for a packet that only spreads.
     */
    void seedGaussian(double cx, double cy, double sigma0,
                      double kx = 0.0, double ky = 0.0) {
        const double denom = 4.0 * sigma0 * sigma0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (onWall(i, j)) { u[i][j] = 0.0; v[i][j] = 0.0; continue; }
                const double x = xAt(j), y = yAt(i);
                const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                const double env = std::exp(-r2 / denom);
                const double phase = kx * x + ky * y;
                u[i][j] = env * std::cos(phase);
                v[i][j] = env * std::sin(phase);
            }
        }
        normalise();
    }

    /**
     * The analytic ground state of the current potential.
     *
     * SquareWell gets sin(pi x/L) sin(pi y/L), which is not merely the
     * continuum ground state but an EXACT eigenvector of the discrete
     * Laplacian with Dirichlet ends - the discrete second difference maps
     * sin(pi j/(N-1)) to (2cos(pi/(N-1)) - 2)/dx^2 times itself. So it is
     * stationary to rounding error, not to O(dx^2), which is what makes it a
     * sharp test case rather than an approximate one.
     *
     * Harmonic gets the Gaussian exp(-omega r^2 / 2), i.e. sigma0 = 1/sqrt(2 omega).
     * That one is only the CONTINUUM eigenstate, so it breathes very slightly
     * on a finite grid.
     *
     * Free has no bound state at all; it gets a broad centred packet, which
     * spreads rather than standing still.
     */
    void seedGroundState() {
        if (pot == Potential::SquareWell) {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    u[i][j] = onWall(i, j)
                        ? 0.0
                        : std::sin(M_PI * i / (rows - 1)) * std::sin(M_PI * j / (cols - 1));
                    v[i][j] = 0.0;
                }
            }
            normalise();
        } else if (pot == Potential::Harmonic) {
            seedGaussian(0.0, 0.0, groundStateSigma());
        } else {
            seedGaussian(0.0, 0.0, 0.1 * cols * dx);
        }
    }

    /**
     * One step: a Cayley transform along x, then along y. Each is separately
     * unitary, so the norm is preserved by construction rather than by luck.
     */
    void step() {
        sweep(true);
        sweep(false);
    }

    /** Probability density |psi|^2 - the field everything downstream plots. */
    const Grid& field() const {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                prob[i][j] = u[i][j] * u[i][j] + v[i][j] * v[i][j];
        return prob;
    }

    /** Total probability. Conserved to rounding error by a unitary scheme. */
    double totalProbability() const {
        double sum = 0.0;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                sum += u[i][j] * u[i][j] + v[i][j] * v[i][j];
        return sum * dx * dx;
    }

private:
    bool onWall(int i, int j) const {
        return i == 0 || j == 0 || i == rows - 1 || j == cols - 1;
    }

    void normalise() {
        const double n = totalProbability();
        if (!(n > 0.0)) return;
        const double s = 1.0 / std::sqrt(n);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) { u[i][j] *= s; v[i][j] *= s; }
    }

    /**
     * One Cayley half-step. Along x the unknowns run across a row; along y,
     * down a column. The coefficients are the same either way, with the
     * direction's own slice of the separable potential, so one function does
     * both rather than two near-identical ones.
     *
     * With s = i dt/2 and k = s / (2 dx^2), the operator (I + s H_dir) is
     *
     *     -k psi_(n-1) + (1 + 2k + s V_n) psi_n - k psi_(n+1)
     *
     * and the right-hand side is the same with s negated. Wall cells solve as
     * identity rows holding psi = 0, the same trick the heat solver uses to
     * carry pinned cells through an implicit solve untouched.
     */
    void sweep(bool alongX) {
        const Complex s(0.0, 0.5 * dt);
        const Complex k = s / (2.0 * dx * dx);
        const int outer = alongX ? rows : cols;
        const int inner = alongX ? cols : rows;

        #ifdef _OPENMP
        #pragma omp parallel
        #endif
        {
        // Per-thread scratch: each grid line is an independent solve, so the
        // sweeps parallelise with no synchronisation, exactly as in the heat
        // solver's ADI.
        std::vector<Complex> a(inner), b(inner), c(inner), d(inner),
                             x(inner), cp(inner), dp(inner);

        #ifdef _OPENMP
        #pragma omp for schedule(runtime)
        #endif
        for (int o = 0; o < outer; ++o) {
            for (int n = 0; n < inner; ++n) {
                const int i = alongX ? o : n;
                const int j = alongX ? n : o;
                if (onWall(i, j)) {
                    a[n] = 0.0; b[n] = 1.0; c[n] = 0.0; d[n] = 0.0;
                    continue;
                }
                const int im = alongX ? i : i - 1, jm = alongX ? j - 1 : j;
                const int ip = alongX ? i : i + 1, jp = alongX ? j + 1 : j;

                const Complex psi(u[i][j], v[i][j]);
                const Complex lo(u[im][jm], v[im][jm]);
                const Complex hi(u[ip][jp], v[ip][jp]);
                const Complex sv = s * (alongX ? vx[j] : vy[i]);

                a[n] = -k;
                b[n] = 1.0 + 2.0 * k + sv;
                c[n] = -k;
                d[n] = (1.0 - 2.0 * k - sv) * psi + k * (lo + hi);
            }
            thomasSolve(a, b, c, d, x, cp, dp);
            for (int n = 0; n < inner; ++n) {
                const int i = alongX ? o : n;
                const int j = alongX ? n : o;
                u[i][j] = x[n].real();
                v[i][j] = x[n].imag();
            }
        }
        } // end parallel region
    }

    int rows, cols;
    double dx, dt;
    Potential pot = Potential::Free;
    double omega = 0.05;
    Grid u, v;               // Re(psi), Im(psi)
    mutable Grid prob;       // |psi|^2, rebuilt on demand so it cannot go stale
    std::vector<double> vx;  // separable potential, V(x,y) = vx[j] + vy[i]
    std::vector<double> vy;
};

#endif // SCHRODINGER_HPP
