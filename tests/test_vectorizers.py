import numpy as np
from microbiome_tool.vectorizers import (
    KmerVectorizer,
    GenomeOceanVectorizer,
    save_vectorizer,
    load_vectorizer,
)

def test_kmer_vectorizer_basic():
    seqs = ["ACGTAC"]
    abundances = [10]
    vec = KmerVectorizer(k=3).transform(seqs, abundances)
    assert isinstance(vec, np.ndarray)
    assert vec.shape[0] > 0
    assert np.isclose(np.linalg.norm(vec), 1.0)


def test_genomeocean_vectorizer_save_load(tmp_path):
    vec = GenomeOceanVectorizer(
        model_name="DOEJGI/GenomeOcean-100M",
        batch_size=4,
        pooling="mean",
        model_max_length=512,
    )
    path = tmp_path / "genomeocean_vec.json"
    save_vectorizer(vec, str(path))

    loaded = load_vectorizer(str(path))
    assert isinstance(loaded, GenomeOceanVectorizer)
    assert loaded.model_name == vec.model_name
    assert loaded.pooling == vec.pooling
    assert loaded.batch_size == vec.batch_size
    assert loaded.model_max_length == vec.model_max_length
