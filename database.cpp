#include "database.hpp"
#include <iostream>
#include <algorithm>

HeatDatabase::HeatDatabase(const std::string& dbName)
    : db(nullptr), insertStmt(nullptr), dbName(dbName), inTransaction(false) {}

HeatDatabase::~HeatDatabase() {
    if (inTransaction) {
        sqlite3_exec(db, "ROLLBACK;", nullptr, nullptr, nullptr);
    }
    if (insertStmt) sqlite3_finalize(insertStmt);
    if (db) sqlite3_close(db);
}

int HeatDatabase::pragmaInt(const char* pragma) {
    sqlite3_stmt* stmt = nullptr;
    int value = -1;
    if (sqlite3_prepare_v2(db, pragma, -1, &stmt, nullptr) == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) value = sqlite3_column_int(stmt, 0);
    }
    sqlite3_finalize(stmt);
    return value;
}

bool HeatDatabase::exec(const char* sql) {
    char* errMsg = nullptr;
    if (sqlite3_exec(db, sql, nullptr, nullptr, &errMsg) != SQLITE_OK) {
        std::cerr << "SQL error [" << sql << "]: "
                  << (errMsg ? errMsg : "unknown") << std::endl;
        sqlite3_free(errMsg);
        return false;
    }
    return true;
}

bool HeatDatabase::init() {
    if (sqlite3_open(dbName.c_str(), &db) != SQLITE_OK) {
        std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
        return false;
    }

    // MUST come first: switching the journal mode writes the database header,
    // and auto_vacuum can no longer be changed once that has happened. This is
    // what lets endRun() hand free pages back to the filesystem instead of
    // letting the file ratchet permanently upward.
    exec("PRAGMA auto_vacuum=INCREMENTAL;");

    // Write-ahead logging with relaxed (but still crash-safe) syncing. The
    // default rollback journal fsyncs on every commit, which dominated runtime
    // when each timestep was its own transaction.
    exec("PRAGMA journal_mode=WAL;");
    exec("PRAGMA synchronous=NORMAL;");
    exec("PRAGMA temp_store=MEMORY;");
    exec("PRAGMA cache_size=-65536;");   // 64 MB page cache

    // One row per saved frame: the whole grid packed row-major as doubles in
    // a single BLOB. Storing one row per cell meant rows*cols B-tree inserts
    // per frame (90,000 at 300x300), which was most of the runtime.
    const char* sql = "CREATE TABLE IF NOT EXISTS Frames ("
                      "step INTEGER PRIMARY KEY, "
                      "rows INTEGER NOT NULL, "
                      "cols INTEGER NOT NULL, "
                      "data BLOB NOT NULL);";
    if (!exec(sql)) return false;

    // Prepared once and reused for every sample of every timestep.
    const char* ins = "INSERT OR REPLACE INTO Frames (step, rows, cols, data) VALUES (?, ?, ?, ?);";
    if (sqlite3_prepare_v2(db, ins, -1, &insertStmt, nullptr) != SQLITE_OK) {
        std::cerr << "Failed to prepare insert: " << sqlite3_errmsg(db) << std::endl;
        return false;
    }
    return true;
}

bool HeatDatabase::beginRun() {
    // Load-bearing only on the shared-database path. Given --run-id, every run
    // opens its own fresh runs/<id>.db and this finds nothing to delete.
    // Without one the binary reuses heat_sim.db across invocations - the CI
    // smoke test and by-hand runs do exactly that - and there the next run
    // would accumulate the previous one's rows and collide with them on the
    // (step, x, y) primary key.
    if (!exec("DELETE FROM Frames;")) return false;
    if (!exec("BEGIN TRANSACTION;")) return false;
    inTransaction = true;
    return true;
}

bool HeatDatabase::saveTimestep(int step, const std::vector<std::vector<double>>& grid) {
    if (!insertStmt) return false;

    const int rows = static_cast<int>(grid.size());
    const int cols = rows > 0 ? static_cast<int>(grid[0].size()) : 0;

    // Flatten into one contiguous buffer, reused across calls.
    frameBuf.resize(static_cast<size_t>(rows) * cols);
    for (int i = 0; i < rows; ++i) {
        std::copy(grid[i].begin(), grid[i].end(),
                  frameBuf.begin() + static_cast<size_t>(i) * cols);
    }

    sqlite3_bind_int(insertStmt, 1, step);
    sqlite3_bind_int(insertStmt, 2, rows);
    sqlite3_bind_int(insertStmt, 3, cols);
    sqlite3_bind_blob(insertStmt, 4, frameBuf.data(),
                      static_cast<int>(frameBuf.size() * sizeof(double)),
                      SQLITE_STATIC);

    bool ok = sqlite3_step(insertStmt) == SQLITE_DONE;
    if (!ok) {
        std::cerr << "Insert failed at step " << step << ": "
                  << sqlite3_errmsg(db) << std::endl;
    }
    sqlite3_reset(insertStmt);
    return ok;
}

bool HeatDatabase::endRun() {
    if (!inTransaction) return true;
    bool ok = exec("COMMIT;");
    inTransaction = false;

    // Release the pages freed by the DELETE in beginRun(). A no-op on a
    // database created before auto_vacuum was enabled.
    exec("PRAGMA incremental_vacuum;");

    // Fallback for such legacy files: once the file is mostly free space, a
    // full VACUUM is the only way to return it. Cheap, because it only runs
    // when the file is already known to be mostly empty.
    int freelist = pragmaInt("PRAGMA freelist_count;");
    int pages    = pragmaInt("PRAGMA page_count;");
    if (pages > 0 && freelist > pages / 2) {
        exec("VACUUM;");
    }

    // Fold the WAL back into the main file so the on-disk size reflects the
    // current run only.
    exec("PRAGMA wal_checkpoint(TRUNCATE);");
    return ok;
}
