import os
import json
from typing import Optional

import numpy as np

from microbiome_tool.utils import read_asv_table, read_asv_table_wide, hash_sequence
from microbiome_tool.vectorizers import (
    KmerVectorizer,
    DNABERT2Vectorizer,
    GenomeOceanVectorizer,
    save_vectorizer,
)
from microbiome_tool.db import (
    DB_FILE,
    init_db,
    insert_asv,
    insert_sample,
    insert_metadata,
)
from microbiome_tool.faiss_index import (
    load_index,
    create_index,
    add_to_index,
    save_index,
)


def load_reference_data(
    asv_table: str,
    metadata: Optional[str] = None,
    *,
    k: int = 5,
    artifacts_dir: str = "artifacts",
    verbose: bool = False,
    vectorizer_type: str = "kmer",
    dnabert2_model_name: Optional[str] = None,
    dnabert2_device: Optional[str] = None,
    genomeocean_model_name: Optional[str] = None,
    genomeocean_model_max_length: Optional[int] = None,
    db_path: Optional[str] = None,
    asv_table_format: str = "long",
):
    """
    Ingest reference ASV data into DuckDB, compute ASV + sample vectors, and update FAISS.

    Steps:
    - Read ASV table with columns: sample_id, asv_id, sequence, abundance
    - Fit KmerVectorizer across all sequences to build a fixed vocabulary
    - Compute ASV-level vectors and persist as .npy (store path in DuckDB)
    - Aggregate per-sample vectors, persist as .npy (store path in DuckDB)
    - Load/create FAISS index, add sample vectors (float32), and persist the index
    - Optionally ingest metadata TSV (per-sample rows) as JSON blobs keyed by sample_id
    """

    def log(msg: str):
        if verbose:
            print(msg, flush=True)

    log(f"[Ingestion] Starting reference ingestion")
    log(f"[Ingestion] ASV table: {asv_table}")
    if metadata:
        log(f"[Ingestion] Metadata file: {metadata}")
    log(f"[Ingestion] Artifacts directory: {artifacts_dir}")
    log(f"[Ingestion] k-mer size: k={k}")

    os.makedirs(artifacts_dir, exist_ok=True)
    asv_vec_dir = os.path.join(artifacts_dir, "asv_vectors")
    sample_vec_dir = os.path.join(artifacts_dir, "sample_vectors")
    mapping_path = os.path.join(artifacts_dir, "index_mapping.json")
    os.makedirs(asv_vec_dir, exist_ok=True)
    os.makedirs(sample_vec_dir, exist_ok=True)

    # 1) Read ASV table
    is_wide = asv_table_format == "wide"
    if is_wide:
        df = read_asv_table_wide(asv_table)
    else:
        df = read_asv_table(asv_table)
    if df.empty:
        raise ValueError("ASV table is empty")

    if is_wide:
        seq_col = df.columns[0]
        n_rows = len(df)
        sample_cols = [c for c in df.columns if c != seq_col]
        n_samples = len(sample_cols)
        log(f"[Step 1/7] Read wide ASV table: {n_rows} sequences, {n_samples} samples")
    else:
        seq_col = "sequence"
        n_rows = len(df)
        n_samples = df["sample_id"].astype(str).nunique()
        log(f"[Step 1/7] Read ASV table: {n_rows} rows, {n_samples} unique samples")

    # 2) Fit vectorizer on all sequences to establish a fixed representation
    vtype = (vectorizer_type or "kmer").lower()
    if vtype == "kmer":
        vectorizer = KmerVectorizer(k=k)
        log(f"[Step 2/7] Fitting KmerVectorizer on {n_rows} sequences (k={k}) ...")
        vectorizer.fit(df[seq_col].astype(str).tolist())
    elif vtype == "dnabert2":
        model_name = dnabert2_model_name or "zhihan1996/DNABERT-2-117M"
        device = None
        if dnabert2_device:
            dnorm = dnabert2_device.lower()
            if dnorm != "auto":
                device = dnorm
        vectorizer = DNABERT2Vectorizer(model_name=model_name, device=device)
        log(
            f"[Step 2/7] Using DNABERT2Vectorizer with model={model_name} "
            f"on {n_rows} sequences ..."
        )
        vectorizer.fit(df[seq_col].astype(str).tolist())
    elif vtype == "genomeocean":
        model_name = genomeocean_model_name or "DOEJGI/GenomeOcean-100M"
        device = None
        if dnabert2_device:
            dnorm = dnabert2_device.lower()
            if dnorm != "auto":
                device = dnorm
        vectorizer = GenomeOceanVectorizer(
            model_name=model_name,
            device=device,
            model_max_length=genomeocean_model_max_length,
        )
        log(
            f"[Step 2/7] Using GenomeOceanVectorizer with model={model_name} "
            f"on {n_rows} sequences ..."
        )
        vectorizer.fit(df[seq_col].astype(str).tolist())
    else:
        raise ValueError(f"Unsupported vectorizer_type: {vectorizer_type}")

    # Persist vectorizer state for future queries (file name kept for backward compat)
    vec_path = os.path.join(artifacts_dir, "kmer_vectorizer.json")
    save_vectorizer(vectorizer, vec_path)
    vocab_size = len(getattr(vectorizer, "vocab_", []) or [])
    log(f"[Step 2/7] Saved vectorizer to {vec_path} (vocab size={vocab_size})")

    # 3) Initialize DB connection
    conn = init_db(db_path=db_path)
    log(f"[Step 3/7] Initialized DuckDB at {db_path or DB_FILE} and ensured tables exist")
    conn.execute("BEGIN TRANSACTION")

    # 4) Optional metadata ingestion: per-sample JSON blobs
    meta_ids = set()
    if metadata:
        import pandas as pd

        log(f"[Step 4/7] Reading metadata TSV ...")
        meta_df = pd.read_csv(metadata, sep="\t")
        for _, row in meta_df.iterrows():
            sample_id = str(row.get("sample_id"))
            if sample_id and sample_id not in meta_ids:
                blob = row.to_dict()
                # Ensure sample_id is explicitly included
                blob["sample_id"] = sample_id
                insert_metadata(conn, sample_id, json.dumps(blob))
                meta_ids.add(sample_id)
        log(f"[Step 4/7] Inserted metadata for {len(meta_ids)} samples")

    # 5) Compute and save ASV vectors; insert ASVs into DB
    # Also collect sample-level inputs for aggregation
    log(f"[Step 5/7] Computing ASV vectors and inserting ASVs into DB ...")
    grouped = {}
    asv_vec_cache = {}
    asv_insert_count = 0

    if not is_wide:
        for _, row in df.iterrows():
            sample_id = str(row["sample_id"])  # type: ignore[index]
            asv_id = str(row["asv_id"])  # type: ignore[index]
            seq = str(row["sequence"])  # type: ignore[index]
            ab = float(row["abundance"])  # type: ignore[index]

            if asv_id not in asv_vec_cache:
                asv_vec = vectorizer.transform([seq], [ab])  # normalized vector
                asv_vec = asv_vec.astype(np.float32)
                asv_vec_path = os.path.join(asv_vec_dir, f"{asv_id}.npy")
                np.save(asv_vec_path, asv_vec)
                asv_vec_cache[asv_id] = asv_vec_path
            else:
                asv_vec_path = asv_vec_cache[asv_id]

            insert_asv(conn, asv_id, seq, ab, sample_id, asv_vec_path)

            bucket = grouped.setdefault(sample_id, {"seqs": [], "abs": []})
            bucket["seqs"].append(seq)
            bucket["abs"].append(ab)
            asv_insert_count += 1
    else:
        seq_col = df.columns[0]
        sample_cols = [c for c in df.columns if c != seq_col]
        for _, row in df.iterrows():
            seq = str(row[seq_col])  # type: ignore[index]
            asv_id = hash_sequence(seq)

            if asv_id not in asv_vec_cache:
                asv_vec = vectorizer.transform([seq], [1.0])
                asv_vec = asv_vec.astype(np.float32)
                asv_vec_path = os.path.join(asv_vec_dir, f"{asv_id}.npy")
                np.save(asv_vec_path, asv_vec)
                asv_vec_cache[asv_id] = asv_vec_path
            else:
                asv_vec_path = asv_vec_cache[asv_id]

            for sample_id in sample_cols:
                val = row[sample_id]
                try:
                    ab = float(val)
                except (TypeError, ValueError):
                    continue
                if ab <= 0:
                    continue

                sample_str = str(sample_id)
                insert_asv(conn, asv_id, seq, ab, sample_str, asv_vec_path)

                bucket = grouped.setdefault(sample_str, {"seqs": [], "abs": []})
                bucket["seqs"].append(seq)
                bucket["abs"].append(ab)
                asv_insert_count += 1

    log(
        f"[Step 5/7] Saved ASV vectors to {asv_vec_dir} "
        f"({len(asv_vec_cache)} unique ASVs, {asv_insert_count} rows)"
    )

    # 6) Aggregate per-sample vectors, persist them, and insert samples into DB
    sample_ids = sorted(grouped.keys())
    sample_vecs = []
    log(f"[Step 6/7] Aggregating sample vectors for {len(sample_ids)} samples ...")
    for sid in sample_ids:
        seqs = grouped[sid]["seqs"]
        abs_ = grouped[sid]["abs"]
        svec = vectorizer.transform(seqs, abs_)
        svec = svec.astype(np.float32)
        spath = os.path.join(sample_vec_dir, f"{sid}.npy")
        np.save(spath, svec)
        # Use sample_id as metadata_id if available, else empty string
        metadata_id = sid if sid in meta_ids or not metadata else ""
        insert_sample(conn, sid, metadata_id, spath)
        sample_vecs.append(svec)
    log(f"[Step 6/7] Saved sample vectors to {sample_vec_dir}")
    conn.execute("COMMIT")

    # 7) Build / update FAISS index and persist
    if not sample_vecs:
        conn.close()
        return {"samples": 0, "asvs": 0}

    dim = sample_vecs[0].shape[0]
    index_path = os.path.join(artifacts_dir, "sample_index.faiss")
    index = load_index(path=index_path)
    mapping = []
    if index is None or getattr(index, "d", dim) != dim:
        # Create a fresh index if none exists or dimension changed
        log(f"[Step 7/7] Creating new FAISS index at {index_path} (dim={dim})")
        index = create_index(dim)
        mapping = []
    else:
        # Load existing mapping if present
        if os.path.exists(mapping_path):
            try:
                with open(mapping_path, "r", encoding="utf-8") as f:
                    mapping = json.load(f)
            except Exception:
                mapping = []
        else:
            mapping = []
        log(f"[Step 7/7] Loaded existing FAISS index and mapping with {len(mapping)} entries")
    mat = np.stack(sample_vecs, axis=0).astype("float32")
    add_to_index(index, mat)
    save_index(index, path=index_path)

    # Update and persist mapping with appended samples in the order added
    mapping.extend(sample_ids)
    with open(mapping_path, "w", encoding="utf-8") as f:
        json.dump(mapping, f)
    log(f"[Step 7/7] Saved FAISS index and updated mapping to {mapping_path} (total={len(mapping)})")

    conn.close()
    summary = {"samples": len(sample_ids), "asvs": asv_insert_count}
    log(f"[Ingestion] Completed: samples={summary['samples']}, asvs={summary['asvs']}")
    return summary

