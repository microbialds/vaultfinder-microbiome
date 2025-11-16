import os
import json
from typing import List, Dict, Optional

import numpy as np

from microbiome_tool.utils import read_asv_table
from microbiome_tool.faiss_index import load_index, search_index
from microbiome_tool.vectorizers import load_vectorizer


def _load_index_mapping(artifacts_dir: str) -> List[str]:
    mapping_path = os.path.join(artifacts_dir, "index_mapping.json")
    if os.path.exists(mapping_path):
        try:
            with open(mapping_path, "r", encoding="utf-8") as f:
                return list(json.load(f))
        except Exception:
            return []
    return []


def compare_query_sample(
    query_table: str,
    top_n: int = 5,
    *,
    artifacts_dir: str = "artifacts",
    db_path: Optional[str] = None,
) -> List[Dict]:
    """
    Transform a query ASV table into a sample-level vector and search the FAISS index.

    Returns a list of dicts with keys: rank, sample_index, sample_id (if available), distance.
    """
    # Load trained vectorizer
    vec_path = os.path.join(artifacts_dir, "kmer_vectorizer.json")
    if not os.path.exists(vec_path):
        # No reference loaded yet; return empty results for permissive behavior
        return []
    vectorizer = load_vectorizer(vec_path)

    # Read and aggregate query ASV table
    df = read_asv_table(query_table)
    if df.empty:
        return []
    seqs = df["sequence"].astype(str).tolist()
    abs_ = df["abundance"].astype(float).tolist()
    qvec = vectorizer.transform(seqs, abs_).astype(np.float32).reshape(1, -1)

    # Load FAISS index
    index_path = os.path.join(artifacts_dir, "sample_index.faiss")
    index = load_index(path=index_path)
    if index is None:
        # No index yet; return empty results
        return []

    distances, indices = search_index(index, qvec, top_n=top_n)
    mapping = _load_index_mapping(artifacts_dir)

    res = []
    for rank, (idx, dist) in enumerate(zip(indices[0], distances[0]), start=1):
        sid = mapping[idx] if mapping and idx < len(mapping) else None
        res.append({
            "rank": rank,
            "sample_index": int(idx),
            "sample_id": sid,
            "distance": float(dist),
        })
    return res
