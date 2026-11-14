"""End-to-end checks of what the real solver binary writes.

tests/api never runs the binary - it seeds databases by hand - so a mismatch
between what database.cpp writes and what server.py reads would pass there.
These tests run ./heat_sim for real and read its output back, both directly
and through the Flask routes the frontend uses.

Needs a compiled ./heat_sim at the repo root (./start.sh or ./run_tests.sh
builds one); the whole module is skipped if it is missing.

    .venv/bin/python -m pytest tests/integration -q
"""
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
import uuid
from array import array

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
BINARY = os.path.join(ROOT, 'heat_sim')
sys.path.insert(0, ROOT)

pytestmark = pytest.mark.skipif(
    not os.access(BINARY, os.X_OK),
    reason='no compiled ./heat_sim - run ./start.sh or ./run_tests.sh first')

import server  # noqa: E402


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """An empty CWD holding a copy of the binary, as server.py expects."""
    shutil.copy2(BINARY, tmp_path / 'heat_sim')
    monkeypatch.chdir(tmp_path)
    return tmp_path


def run_solver(*args):
    """Run the binary with a fresh run id; return (run_id, frames, json_out).

    frames is [(step, rows, cols, grid)] in step order, unpacked here
    independently of server.py so the two readers can be compared.
    """
    run_id = uuid.uuid4().hex
    proc = subprocess.run(['./heat_sim', *map(str, args), f'--run-id={run_id}'],
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr or proc.stdout

    conn = sqlite3.connect(f'runs/{run_id}.db')
    rows = conn.execute('SELECT step, rows, cols, data FROM Frames ORDER BY step').fetchall()
    conn.close()

    frames = []
    for step, r, c, blob in rows:
        assert len(blob) == r * c * 8, f'step {step}: {len(blob)} bytes for {r}x{c}'
        flat = array('d')
        flat.frombytes(blob)
        frames.append((step, r, c, [flat[i * c:(i + 1) * c].tolist() for i in range(r)]))

    with open(f'runs/{run_id}.json') as f:
        out = json.load(f)
    return run_id, frames, out


# Non-square on purpose: a rows/cols swap anywhere in the pipeline breaks it.
ROWS, COLS = 23, 37


@pytest.mark.parametrize('mode,F', [('fdm', 0.2), ('be', 5.0), ('cn', 5.0)])
def test_every_frame_has_the_right_shape(workdir, mode, F):
    _, frames, _ = run_solver(ROWS, COLS, 100, 0, 0, 0, mode, 0.01, F)
    assert len(frames) >= 2
    for step, r, c, grid in frames:
        assert (r, c) == (ROWS, COLS), f'step {step}'
        assert len(grid) == ROWS and all(len(row) == COLS for row in grid)


@pytest.mark.parametrize('mode,F', [('fdm', 0.2), ('be', 5.0), ('cn', 5.0)])
def test_last_frame_matches_the_final_json_grid(workdir, mode, F):
    # The JSON is written straight from the solver's grid, independently of
    # the database path, so the last stored frame must be the same grid. The
    # JSON prints 4 decimals (exportToJSON's setprecision(4)); the database
    # keeps full doubles, so each stored value must round to the JSON one.
    _, frames, out = run_solver(ROWS, COLS, 100, 0, 0, 0, mode, 0.01, F)
    step, _, _, grid = frames[-1]
    assert step == out['step']
    for i, (db_row, js_row) in enumerate(zip(grid, out['data'])):
        for j, (db, js) in enumerate(zip(db_row, js_row)):
            assert abs(db - js) <= 5e-5 + 1e-12, f'cell ({i},{j}): {db} vs {js}'
    # Not a trivial grid: the interior has actually been solved.
    assert 0.0 < grid[ROWS // 2][COLS // 2] < 100.0


def test_steps_start_at_zero_and_increase(workdir):
    _, frames, out = run_solver(ROWS, COLS, 100, 0, 0, 0, 'cn', 0.01, 5)
    steps = [f[0] for f in frames]
    assert steps[0] == 0
    assert steps == sorted(set(steps))
    assert steps[-1] == out['step']


def test_boundaries_land_on_the_right_edges_in_every_frame(workdir):
    # Distinct temperature per edge: top=100, bottom=-20, left=35, right=7.
    # Corners belong to whichever edge main.cpp writes last, so skip them.
    _, frames, _ = run_solver(ROWS, COLS, 100, -20, 35, 7, 'cn', 0.01, 5)
    for step, _, _, g in frames:
        assert all(v == 100.0 for v in g[0][1:-1]), f'top, step {step}'
        assert all(v == -20.0 for v in g[-1][1:-1]), f'bottom, step {step}'
        assert all(row[0] == 35.0 for row in g[1:-1]), f'left, step {step}'
        assert all(row[-1] == 7.0 for row in g[1:-1]), f'right, step {step}'


@pytest.mark.parametrize('mode,F', [('fdm', 0.2), ('cn', 5.0)])
def test_a_heat_source_stays_pinned_in_every_frame(workdir, mode, F):
    # Trailing (row, col, temp) triple; off-centre so a transpose misses it.
    r, c = 6, 25
    _, frames, _ = run_solver(ROWS, COLS, 0, 0, 0, 0, mode, 0.01, F, r, c, 500)
    for step, _, _, g in frames:
        assert g[r][c] == 500.0, f'step {step}'
    # ...and heat actually spreads from it, so the frames are not all identical.
    assert frames[-1][3][r][c + 1] > frames[0][3][r][c + 1]


def test_server_routes_serve_exactly_what_the_binary_stored(workdir, monkeypatch):
    """POST /run through Flask with the real binary, then read it back.

    This is the path the browser takes: the run id from the POST response,
    /frames for the slider, /run?time= for each frame.
    """
    monkeypatch.setattr(server.limiter, 'enabled', False)
    client = server.app.test_client()

    res = client.post('/run', json={'rows': ROWS, 'cols': COLS, 'mode': 'cn',
                                    'top': 100, 'bottom': 0, 'left': 0, 'right': 0})
    assert res.status_code == 200, res.get_data(as_text=True)
    body = res.get_json()
    run_id = body.get('run_id') or body['data']['run_id']

    steps = client.get(f'/frames?run_id={run_id}').get_json()['steps']
    assert steps and steps[0] == 0

    conn = sqlite3.connect(f'runs/{run_id}.db')
    stored = dict(((s, (r, c, d)) for s, r, c, d in
                   conn.execute('SELECT step, rows, cols, data FROM Frames')))
    conn.close()
    assert sorted(stored) == steps

    # First, a middle one, and the last - each must match the raw bytes.
    for step in {steps[0], steps[len(steps) // 2], steps[-1]}:
        frame = client.get(f'/run?time={step}&run_id={run_id}').get_json()
        r, c, blob = stored[step]
        flat = array('d', [v for row in frame['data'] for v in row])
        assert (frame['step'], frame['rows'], frame['cols']) == (step, r, c)
        assert flat.tobytes() == blob

    # A time between two stored frames snaps down to the earlier one.
    if len(steps) >= 2 and steps[1] - steps[0] > 1:
        mid = client.get(f'/run?time={steps[1] - 1}&run_id={run_id}').get_json()
        assert mid['step'] == steps[0]


# ---- Regression guards on the speed-up itself ----------------------------
# Generous bounds: they exist to catch a return to per-cell storage (5 s and
# 175 MB before the change), not to benchmark the machine.

def test_explicit_100x100_finishes_quickly(workdir):
    start = time.perf_counter()
    run_solver(100, 100, 100, 0, 0, 0, 'fdm', 0.01, 0.2)
    elapsed = time.perf_counter() - start
    assert elapsed < 2.5, f'took {elapsed:.2f}s (was ~5s with one row per cell)'


def test_database_stays_compact(workdir):
    run_id, frames, _ = run_solver(100, 100, 100, 0, 0, 0, 'fdm', 0.01, 0.2)
    size = os.path.getsize(f'runs/{run_id}.db')
    payload = len(frames) * 100 * 100 * 8
    # One row per frame costs a little over the raw doubles; one row per cell
    # cost about 2.7x them.
    assert size < 1.5 * payload, f'{size / 1e6:.1f} MB for {payload / 1e6:.1f} MB of data'
