import duckdb

DB_FILE = "reference.duckdb"

def init_db(db_path=None):
    """Initialize (or connect to) DuckDB and ensure required tables exist."""
    path = db_path or DB_FILE
    conn = duckdb.connect(path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS samples (
            sample_id STRING,
            metadata_id STRING,
            sample_vector_path STRING
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS asvs (
            asv_id STRING,
            sequence STRING,
            abundance FLOAT,
            sample_id STRING,
            vector_path STRING
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            metadata_id STRING,
            json_blob JSON
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS queries (
            query_id STRING,
            timestamp TIMESTAMP,
            input_file STRING,
            results_file STRING
        )
    """)
    return conn

def insert_sample(conn, sample_id, metadata_id, sample_vector_path):
    conn.execute(
        "INSERT INTO samples VALUES (?, ?, ?)",
        [sample_id, metadata_id, sample_vector_path]
    )

def insert_asv(conn, asv_id, seq, abundance, sample_id, vector_path):
    conn.execute(
        "INSERT INTO asvs VALUES (?, ?, ?, ?, ?)",
        [asv_id, seq, abundance, sample_id, vector_path]
    )

def insert_metadata(conn, metadata_id, json_blob):
    conn.execute(
        "INSERT INTO metadata VALUES (?, ?)",
        [metadata_id, json_blob]
    )
