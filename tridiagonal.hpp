#ifndef TRIDIAGONAL_HPP
#define TRIDIAGONAL_HPP

#include <vector>

/**
 * Thomas algorithm: O(n) solve of a tridiagonal system.
 *
 * Extracted from HeatSimulation so the Schrodinger solver can use the same
 * code. The scalar type is the only thing that differs between them: the heat
 * ADI sweeps solve real systems, while Crank-Nicolson on a complex wave
 * function solves complex ones. The algorithm itself only ever uses +, -, *
 * and /, so one template covers both.
 *
 * Only this scalar solve is shared. The surrounding ADI sweeps are not: the
 * heat coefficients come from harmonic-mean face conductivities and pinned
 * boundary rows, the Schrodinger ones from a complex Cayley transform with a
 * potential on the diagonal. They have no common structure worth abstracting.
 *
 * No pivoting, which is safe only for a diagonally dominant system - true for
 * both callers here. The heat sweeps give b = 1 + rc(aW+aE) against
 * |a| + |c| = rc(aW+aE); the Schrodinger sweeps give |b| >= sqrt(1 + 4k^2)
 * against |a| + |c| = 2k for any potential V >= 0. A potential allowed to go
 * negative would need that checked again.
 *
 * cp and dp are caller-supplied scratch of the same length, so a sweep can
 * allocate once per thread rather than once per grid line.
 */
template <typename T>
void thomasSolve(const std::vector<T>& a, const std::vector<T>& b,
                 const std::vector<T>& c, const std::vector<T>& d,
                 std::vector<T>& x,
                 std::vector<T>& cp, std::vector<T>& dp) {
    const int n = static_cast<int>(d.size());
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        const T m = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
    }
    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = dp[i] - cp[i] * x[i + 1];
}

#endif // TRIDIAGONAL_HPP
