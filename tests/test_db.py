from microbiome_tool import db

def test_init_db():
    conn = db.init_db()
    tables = conn.execute("SHOW TABLES").fetchall()
    expected = {"samples", "asvs", "metadata", "queries"}
    got = {t[0] for t in tables}
    assert expected.issubset(got)


def test_init_db_custom_path(tmp_path):
    custom_db = tmp_path / "custom.duckdb"
    conn = db.init_db(db_path=str(custom_db))
    tables = conn.execute("SHOW TABLES").fetchall()
    expected = {"samples", "asvs", "metadata", "queries"}
    got = {t[0] for t in tables}
    assert expected.issubset(got)
    assert custom_db.exists()
