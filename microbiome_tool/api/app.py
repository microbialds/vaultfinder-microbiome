from fastapi import FastAPI, UploadFile
import pandas as pd
import numpy as np
import os
import json
from typing import Optional
from microbiome_tool.vectorizers import load_vectorizer
from microbiome_tool.faiss_index import load_index, search_index

app = FastAPI(title="Microbiome Comparison API")

@app.get("/")
def root():
    return {"message": "Microbiome Comparison API is running"}

@app.post("/compare")
async def compare_sample(
    file: UploadFile,
    top_n: int = 5,
    artifacts_dir: str = "artifacts",
    db_path: Optional[str] = None,
):
    """Compare a query ASV table (TSV) against the FAISS index using saved vectorizer."""
    content = await file.read()
    df = pd.read_csv(pd.io.common.StringIO(content.decode()), sep="\t")
    if df.empty:
        return {"results": []}

    vec_path = os.path.join(artifacts_dir, "kmer_vectorizer.json")
    if not os.path.exists(vec_path):
        return {"error": "No saved vectorizer found. Please load reference data first."}
    vectorizer = load_vectorizer(vec_path)

    seqs = df["sequence"].astype(str).tolist()
    abs_ = df["abundance"].astype(float).tolist()
    query_vec = vectorizer.transform(seqs, abs_).astype("float32").reshape(1, -1)

    index_path = os.path.join(artifacts_dir, "sample_index.faiss")
    index = load_index(path=index_path)
    if index is None:
        return {"error": "No FAISS index found. Please load reference data first."}
    distances, indices = search_index(index, query_vec, top_n=top_n)

    # Load mapping of index positions to sample_ids
    mapping_path = os.path.join(artifacts_dir, "index_mapping.json")
    mapping = []
    if os.path.exists(mapping_path):
        try:
            with open(mapping_path, "r", encoding="utf-8") as f:
                mapping = json.load(f)
        except Exception:
            mapping = []

    results = []
    for i, (idx, dist) in enumerate(zip(indices[0], distances[0])):
        sid = mapping[idx] if mapping and idx < len(mapping) else None
        results.append({
            "rank": i + 1,
            "sample_index": int(idx),
            "sample_id": sid,
            "distance": float(dist),
        })
    return {"results": results}
