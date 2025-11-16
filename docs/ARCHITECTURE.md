# Architecture

This document outlines the core components and data flow of the microbiome comparison tool.

## Components

- `microbiome_tool/vectorizers.py`
  - `Vectorizer` base class.
  - `KmerVectorizer` with `fit()` to build a sorted k-mer vocabulary and `transform()` to produce L2-normalized vectors.
  - `DNABERT2Vectorizer` and `GenomeOceanVectorizer` for transformer-based embeddings.
  - JSON helpers to save/load vectorizer parameters (k-mer vocabulary or transformer configuration), used by the ingestion and comparison pipelines.

- `microbiome_tool/db.py`
  - Initializes DuckDB (`reference.duckdb`)
  - Tables:
    - `samples(sample_id, metadata_id, sample_vector_path)`
    - `asvs(asv_id, sequence, abundance, sample_id, vector_path)`
    - `metadata(metadata_id, json_blob)`
    - `queries(query_id, timestamp, input_file, results_file)`

- `microbiome_tool/faiss_index.py`
  - `create_index(dim)`, `add_to_index(vectors)`, `search_index(query_vec, top_n)`, `save_index(path)`, `load_index(path)`

- `microbiome_tool/pipeline/ingestion.py`
  - `load_reference_data(asv_table, metadata=None, *, k=5, artifacts_dir="artifacts", vectorizer_type="kmer", ...)`.
  - Reads ASV TSV, selects and fits a vectorizer (`KmerVectorizer`, `DNABERT2Vectorizer`, or `GenomeOceanVectorizer`), writes ASV and sample vectors to `.npy` files, inserts into DuckDB, builds/updates FAISS index, and persists index mapping.

- `microbiome_tool/pipeline/comparison.py`
  - `compare_query_sample(query_table, top_n=5, *, artifacts_dir="artifacts")`
  - Loads vectorizer + FAISS, computes query vector, searches, returns ranked matches with distances and sample IDs via mapping.

- `microbiome_tool/api/app.py`
  - FastAPI endpoints: `/` and `/compare`
  - `/compare` loads vectorizer + FAISS from `artifacts/` and returns ranked matches.

- `cli/`
  - `load_reference` (builds artifacts)
  - `compare_samples` (searches FAISS)

## Data flow: ingestion

1. Read ASV table (`sample_id, asv_id, sequence, abundance`).
2. Fit `KmerVectorizer(k)` on all sequences to build a deterministic vocabulary.
3. For each ASV row, compute and save a normalized ASV-level vector (`.npy`), record in DuckDB `asvs`.
4. For each sample, aggregate sequences+abundances into a sample-level vector, save (`.npy`), record in DuckDB `samples`.
5. Build or update FAISS index with sample vectors, persist to `artifacts/sample_index.faiss`.
6. Persist ordered mapping from FAISS row index to `sample_id` at `artifacts/index_mapping.json`.
7. Optionally ingest metadata TSV rows into DuckDB `metadata` as JSON blobs keyed by `sample_id`.

## Data flow: comparison

1. Load `kmer_vectorizer.json` and `sample_index.faiss` from `artifacts/`.
2. Read query ASV TSV, compute query vector using the saved vocabulary.
3. Search FAISS for `top_n` nearest neighbors.
4. Map FAISS row indices to `sample_id` via `index_mapping.json`.
5. Return ranked results with distances.

## Reproducibility

- The sorted k-mer vocabulary ensures fixed vector dimensionality and ordering across runs.
- Vector normalization is applied for consistent distance computations.
- FAISS index and vectorizer state are persisted and reused between runs.

## Extensibility

- To add a new vectorizer (e.g., DNA-BERT):
  - Implement a subclass of `Vectorizer` with `fit()` and `transform()`.
  - Provide persistence (save/load) helpers.
  - Swap it into the ingestion and comparison workflows.

## File locations

- Artifacts directory: `artifacts/`
- DuckDB database: `reference.duckdb`
- CLI entrypoints: `load-reference`, `compare-samples`
- API: `uvicorn microbiome_tool.api.app:app --reload`
