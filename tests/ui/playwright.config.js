// @ts-check
const { defineConfig } = require('@playwright/test');

// The suite is driven against a server that run_ui_tests.sh (or CI) starts
// separately - this config never launches the app itself, since building and
// starting it needs the C++ binary compiled first.
module.exports = defineConfig({
  testDir: './tests',
  // Not the default 5s: a run has to launch the solver and come back.
  expect: { timeout: 15_000 },
  // Conservative now, rather than required. This was set when every run shared
  // one database, so a GET after this suite's own POST could land after a
  // different worker's run had already cleared the table. Per-run files (see
  // RUNS_DIR in server.py) removed that race. What is left is weaker: the rate
  // limiter counts per IP across the whole suite, and SIM_LOCK serialises the
  // solver processes anyway, so parallel workers would buy little. Check both
  // before raising this.
  workers: 1,
  retries: process.env.CI ? 1 : 0,
  reporter: process.env.CI ? [['list'], ['html', { open: 'never' }]] : 'list',
  use: {
    baseURL: process.env.UI_BASE_URL || 'http://127.0.0.1:5000',
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
  },
});
