from microbiome_tool.pipeline import ingestion


def test_wide_ingestion_builds_artifacts_and_counts_asvs(tmp_path):
    artifacts = tmp_path / "artifacts_wide"

    # Create a small wide-format ASV table locally
    wide_path = tmp_path / "asv_table_wide.tsv"
    wide_path.write_text(
        "sequence\tS1\tS2\n"  # header
        "ACGTACGTAC\t10\t0\n"  # present in S1 only
        "TGCACTGACT\t0\t5\n",  # present in S2 only
        encoding="utf-8",
    )

    # Ingest the small wide-format ASV table
    res = ingestion.load_reference_data(
        asv_table=str(wide_path),
        metadata="data/meta.tsv",
        k=3,
        artifacts_dir=str(artifacts),
        asv_table_format="wide",
    )

    # Basic shape of the summary
    assert "samples" in res
    assert "asvs" in res
    assert res["samples"] > 0
    assert res["asvs"] > 0

    # Check that core artifacts were created
    asv_vec_dir = artifacts / "asv_vectors"
    sample_vec_dir = artifacts / "sample_vectors"
    index_file = artifacts / "sample_index.faiss"
    mapping_file = artifacts / "index_mapping.json"

    assert asv_vec_dir.exists() and any(asv_vec_dir.iterdir())
    assert sample_vec_dir.exists() and any(sample_vec_dir.iterdir())
    assert index_file.exists()
    assert mapping_file.exists()
