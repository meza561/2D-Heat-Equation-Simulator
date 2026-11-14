# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
./start.sh          # build heat_sim, set up the venv, launch Flask on :5000, open browser
./run_tests.sh       # build + run heat_tests (54 physics assertions), exits non-zero on failure
./run_study.sh        # build + run heat_study (convergence/stability/scaling), render figures/ via make_figures.py
./run_ui_tests.sh     # build heat_sim, start a server, run the Playwright suite headless, stop the server
venv/bin/python -m pytest tests/api -q   # Flask API tests (59): no server, no binary, ~0.1s
```

The scripts create `.venv` (they accept a pre-existing `venv` too, but `.venv` is
what they make) — **this checkout has `venv/`**, which is why the command above
uses that path. They install `requirements.txt` only. The other two
requirement files are opt-in, one-time installs kept out of the runtime set so
the Docker image ships neither: `venv/bin/pip install -r requirements-test.txt`
for pytest, which `tests/api` needs and no script installs for you, and
`requirements-study.txt` for matplotlib, which `run_study.sh` prompts for when
the figures step is skipped. CI installs the test one explicitly.

All four shell scripts prefer CMake (`build/` dir) but fall back to compiling the relevant `.cpp` directly with `clang++`/`g++` if no CMake is found (they search `PATH` then CLion's bundled copy). After a CMake build, the binary is copied over the top-level `./heat_sim` / `./heat_study` (atomically, via rename) because `server.py` and ad-hoc runs invoke it from the project root — and it's re-signed with `codesign --force --sign -` since macOS SIGKILLs a binary whose ad-hoc signature doesn't match its inode after an in-place copy.

`run_ui_tests.sh` binds port 5000 like `start.sh` does, and refuses to start if something else holds it; `PORT=5577 ./run_ui_tests.sh` runs against another port instead of fighting a dev instance (`server.py` reads `PORT`, defaulting to 5000).

There is no separate lint step. The C++ test binary isn't structured for filtering by name, so to run a single physics assertion it's easiest to isolate it in `tests.cpp` and re-run `./run_tests.sh`; the other two suites filter natively (`pytest -k <substring>`, `npx playwright test <file> -g <title>` from `tests/ui/`).

To invoke the simulator directly instead of through the API:
```bash
./heat_sim ROWS COLS top bottom left right [mode] [alpha] [F] [row col temp]...
# mode: fdm (explicit, default) | be (backward Euler) | cn (Crank-Nicolson) | pde (analytical) | fisher | gray-scott | wave | schrodinger
# optional flags (anywhere in argv): --insulate=tblr   --material=r0,c0,r1,c1,alpha (repeatable)
```
Fisher-KPP/Gray-Scott/wave/Schrodinger take their own trailing params instead of boundary temps/sources — see `main.cpp`'s argv parsing for the exact positions.

## Architecture

**Three-tier, one run per file.** `server.py` (Flask) shells out to the compiled `./heat_sim` binary per request (`subprocess.run`, 300s timeout) rather than linking it as a library. Every run gets its own `runs/<run_id>.db` and `runs/<run_id>.json`, named by a uuid4 the `/run` handler generates and passes to the binary as `--run-id=`. That replaced an earlier design where every run wrote one shared `heat_sim.db` and `latest_heatmap.json`, which is why the binary still accepts no `--run-id` and falls back to that shared pair (the CI smoke test and by-hand invocations rely on it).

`SIM_LOCK` survives that change but no longer protects anything: runs cannot corrupt each other now, so it is purely a cap on how many CPU-bound solver processes start at once. A `run_id` reaches the filesystem, so `_valid_run_id` rejects anything that is not 32 lowercase hex of a version-4 uuid, and `/frames` and `GET /run?time=` both **require** one — 400 otherwise. `sweep_old_runs()` deletes run files past `RUN_RETENTION_HOURS` (default 24) at import and before each run, so `runs/` cannot grow without bound.

**Coordinate flip happens once, at the API boundary.** The UI and API use bottom-left-origin (x right, y up); the solver uses `[row][col]` with row 0 at the top. `server.py`'s `/run` handler is the only place this conversion occurs (for point sources and material rectangles) — `main.cpp` and `simulation.hpp` are unaware of the UI's coordinate convention.

**`simulation.hpp` is the shared solver core**, linked into `heat_sim` (the app), `heat_study` (validation/scaling), and `heat_tests` (regression suite) — all three targets in `CMakeLists.txt` compile it directly rather than linking a shared library, so a change there affects all three binaries and should be re-validated with `./run_tests.sh` and usually `./run_study.sh`. `reaction.hpp` (Fisher-KPP / Gray-Scott), `wave.hpp` (2D wave equation) and `schrodinger.hpp` (2D time-dependent Schrodinger) are parallel, separate models sharing only the five-point Laplacian idea, not code, with the heat solver. The one genuine code share is `tridiagonal.hpp`: the Thomas solver, templated on the scalar type so the heat ADI sweeps (real) and the Schrodinger Cayley sweeps (complex) run the same routine. The sweeps around it are not shared and should not be — heat's coefficients come from harmonic-mean face conductivities and pinned rows, Schrodinger's from a complex Cayley transform with a potential on the diagonal. `wave.hpp` is the one that is second order in time, so it carries two displacement grids and starts with a half-coefficient Taylor step rather than a leapfrog one — that factor is what the observed-order assertion in `tests.cpp` protects. `schrodinger.hpp` is the one whose state is complex; it stores psi as two real grids so the export path, database schema and frontend keep taking real matrices, and exports `|psi|^2`.

**Time integration.** Explicit (forward Euler) is conditionally stable (F = alpha·dt/dx² ≤ 0.25); Backward Euler and Crank-Nicolson use ADI (alternating-direction-implicit) splitting — Douglas-Rachford and Peaceman-Rachford respectively — reducing each 2D step to one tridiagonal (Thomas algorithm) solve per grid line, and are unconditionally stable so `F` becomes a free UI parameter for them. `dt` is always derived from `F` and `alpha`, never fixed, because a fixed dt at the default alpha under-resolved the approach to steady state (see `main.cpp` comments). Dirichlet boundaries and pinned point sources enter the implicit systems as identity rows.

**Heterogeneous media** use divergence-form diffusion (`∇·(α∇u)`, not `α∇²u`) with harmonic-mean conductivity at cell faces — required to get flux right across a material discontinuity; this is load-bearing for the composite-slab regression test.

**Persistence contract:** `database.cpp`/`.hpp` wraps SQLite, one database per run. `beginRun()` opens with `DELETE FROM HeatMap`, which is a no-op for a per-run file but still load-bearing on the shared-database fallback: without `--run-id` the binary reuses `heat_sim.db` across invocations, where the next run would otherwise accumulate the previous one's rows and collide on the `(step, x, y)` primary key. Timesteps are saved every `SAVE_INTERVAL` steps, and the save interval widens automatically on large grids to bound file size under a fixed sample budget — the frontend's `/frames?run_id=` endpoint lists the actually-stored step numbers so the time slider addresses real frames instead of guessing from the interval, and `/run?time=X&run_id=` snaps down to the nearest stored frame. The frontend keeps the `run_id` from the POST response and sends it with every subsequent frame request; history entries store it too.

**Parallelism:** step loops use `#pragma omp` with `schedule(runtime)` (default `static`; compare via `OMP_SCHEDULE=guided ./run_study.sh`). OpenMP is optional and compiles away cleanly when absent (Apple's clang doesn't ship it — `CMakeLists.txt` has custom logic to locate Homebrew's `libomp` on macOS since CMake can't find it unaided). Parallel and serial runs must stay bit-identical — rows of the explicit update and each ADI grid line are independent, so no synchronization is needed; preserve this property when touching the step loops.

**Validation is data-driven, not just pass/fail:** `heat_study` writes `study_results.json` + per-study CSVs, which `server.py`'s `/study` endpoint serves and `web/validation.js` renders in the app's Validation tab. `make_figures.py` renders the same data into `figures/*.png`/`.pdf` for the README. If you change solver numerics, re-run `./run_study.sh` and check whether the measured orders/thresholds in `README.md`'s tables still hold.

## Tests

Three suites, three CI jobs (`.github/workflows/ci.yml`), each covering a different layer:

- **`tests.cpp` → `heat_tests`** (54 checks) — the physics. Every assertion is a value theory fixes independently of the code (exact analytical temperatures, convergence orders, the F = 1/4 threshold, conservation laws, the wave CFL limit, Schrodinger norm conservation and the free-packet spreading law), so it's the suite to re-run after any change to `simulation.hpp`, `reaction.hpp`, `wave.hpp`, `schrodinger.hpp` or `tridiagonal.hpp`.
- **`tests/api/` (pytest, 59 tests)** — the Flask layer in isolation: clamping, `/frames`, `GET /run?time=`'s snap-to-nearest-earlier-frame, `POST /run`'s error branches with `subprocess.run` mocked, per-run isolation and `run_id` validation (including a path-traversal attempt), the retention sweep, the rate limiter's JSON 429, and the deployment routes and security headers. Needs no server and no compiled binary. The fixture `chdir`s into a tmp dir, because `server.py` resolves `runs/` relative to the working directory — keep it that way, or tests start writing into the real `runs/`. It also disables the limiter, whose counter is process-global and would otherwise be consumed across unrelated tests. Because the handler generates the `run_id` itself, the mocked `subprocess.run` parses `--run-id=` out of argv and writes the JSON where the real binary would.
- **`tests/ui/` (Playwright, Node, 19 tests)** — the frontend end to end against a real server: every simulation mode (including one per Schrodinger potential), the timeline indexing real stored frames by `run_id`, the legend gradient matching what the canvas paints, hover mapping on a non-square grid, history round-tripping, console errors, and narrow-viewport overflow.

Two constraints in the UI suite that are easy to undo by accident:

- **`workers: 1`, though the original reason is gone.** It was set when every run shared one database, so a GET after one test's POST could land after another test's run had cleared the table. Per-run files removed that race. What is left is weaker: the rate limiter counts per IP across the whole suite, and `SIM_LOCK` serialises the solver anyway so parallel workers buy little. Treat it as conservative rather than load-bearing — but verify against the limiter and lock queueing before raising it.
- **Assertions compare the app's own outputs to each other, not to a reimplementation.** The legend test samples the canvas bitmap and the rendered CSS gradient and checks they agree; duplicating `heatColor()` in the test would just be a second copy free to drift the same way the old hardcoded gradient did.

`study_results.json` is a gitignored `./run_study.sh` output that CI never has, so `responsive.spec.js` serves `tests/ui/fixtures/study_results.json` through `page.route` to get the Validation tab's charts rendered at all (measuring before they load is how a 76px overflow went unnoticed). That fixture needs its `.gitignore` negation to stay in place — the root `study_results.json` pattern is unanchored and otherwise matches it too. `console-errors.spec.js` deliberately doesn't mock it: it waits for whichever state the real server produces, so the "not generated yet" path stays covered.

## Deployment

Containerised: `Dockerfile` builds the solver and runs `gunicorn wsgi:app`;
`fly.toml` is a worked example, though the live site runs on Render behind
Cloudflare. **One worker, deliberately** — `wsgi.py` explains why: the rate
limiter counts in process memory, and N workers means N concurrent solvers.
Environment variables are documented in `.env.example`; see `README.md`'s
"Deployment" section for the commands.

## Project structure

See `README.md`'s "Project Structure" section and its "API Endpoints" section for the `/run` payload shape — both are current and detailed; don't duplicate them here.
