import os
import json
from microbiome_tool.pipeline import ingestion


def test_ingestion_runs_and_persists(tmp_path):
    artifacts = tmp_path / "artifacts"
    res = ingestion.load_reference_data(
        asv_table="data/asv_table.tsv",
        metadata="data/meta.tsv",
        k=3,
        artifacts_dir=str(artifacts),
    )
    assert res["samples"] == 3
    assert res["asvs"] == 9

    # Check artifacts exist
    assert (artifacts / "kmer_vectorizer.json").exists()
    assert (artifacts / "sample_index.faiss").exists()
    assert (artifacts / "index_mapping.json").exists()
    # Mapping should have 3 sample IDs
    mapping = json.loads((artifacts / "index_mapping.json").read_text())
    assert len(mapping) == 3
