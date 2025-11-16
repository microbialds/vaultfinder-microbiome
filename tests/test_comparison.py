from microbiome_tool.pipeline import comparison, ingestion


def test_compare_returns_empty_when_no_artifacts(tmp_path):
    # Use isolated artifacts dir with no vectorizer/index
    results = comparison.compare_query_sample(
        query_table="data/asv_table.tsv",
        top_n=2,
        artifacts_dir=str(tmp_path / "artifacts"),
    )
    assert results == []


def test_compare_after_ingestion(tmp_path):
    artifacts = tmp_path / "artifacts"
    # Ingest reference first
    ingestion.load_reference_data(
        asv_table="data/asv_table.tsv",
        metadata="data/meta.tsv",
        k=3,
        artifacts_dir=str(artifacts),
    )
    # Then compare
    results = comparison.compare_query_sample(
        query_table="data/asv_table.tsv",
        top_n=2,
        artifacts_dir=str(artifacts),
    )
    assert isinstance(results, list)
    assert len(results) == 2
    # Ensure required fields are present
    r0 = results[0]
    assert set(["rank", "sample_index", "distance"]).issubset(r0.keys())
