# Usage Guide

This guide explains how to ingest reference ASV tables, build vectors and a FAISS index, and run queries via CLI, API, or Python.

## Requirements

- Python >= 3.10
- Packages are listed in `requirements.txt` and installed by `pip install -e .`
- CPU-only FAISS is used (`faiss-cpu`).

## Data formats

- ASV table (TSV), either:
  - Long format (default) with the following columns:
    - `sample_id`
    - `asv_id`
    - `sequence`
    - `abundance`
  - Wide format, where the first column is the DNA sequence and each subsequent column is a sample ID with the abundance of that sequence in the sample.
- Metadata (TSV, optional): must contain `sample_id` and any other columns you wish to store.

Example long-format ASV table: `data/asv_table.tsv`

```
sample_id	asv_id	sequence	abundance
S1	ASV1	ACGTACGTAC	50
...
```

Example wide-format ASV table (sequence x samples):

```
sequence	S1	S2	S3
ACGTACGTAC	50	40	25
TGCACTGACT	30	25	0
...
```

Wide-format tables are currently supported for reference ingestion; query ASV tables should use the long format.

## Quickstart (CLI)

1) Load reference data and build artifacts (vectorizer, vectors, FAISS, mapping):

```
load-reference --asv_table data/asv_table.tsv --metadata data/meta.tsv
```

To ingest a wide-format ASV table (sequence x samples):

```
load-reference --asv_table data/asv_table_wide.tsv --metadata data/meta.tsv --asv_table_format wide
```

2) Run a comparison of a query ASV table against the FAISS index:

```
compare-samples --query data/asv_table.tsv --top_n 5
```

By default, these commands use a DuckDB database called `reference.duckdb` at the
repo root. To write to a different database file when ingesting, pass
`--db_path my_reference.duckdb` to `load-reference`. When comparing, you can
also pass `--artifacts_dir` and `--db_path` if you manage multiple references
with separate artifacts and databases.

To build the reference using DNABERT2 embeddings instead of k-mers, use the
`dnabert2` vectorizer type and (optionally) provide a model name:

```
load-reference \
  --asv_table data/asv_table.tsv \
  --metadata data/meta.tsv \
  --vectorizer_type dnabert2 \
  --dnabert2_model_name zhihan1996/DNABERT-2-117M
```

To build the reference using GenomeOcean-100M (Mistral-based) embeddings, use
the `genomeocean` vectorizer type and (optionally) provide a model name and
maximum token length. The default model is `DOEJGI/GenomeOcean-100M`:

```
load-reference \
  --asv_table data/asv_table.tsv \
  --metadata data/meta.tsv \
  --vectorizer_type genomeocean \
  --genomeocean_model_name DOEJGI/GenomeOcean-100M \
  --genomeocean_model_max_length 1024 \
  --device cuda
```

Output includes the FAISS row index and the resolved `sample_id` when available.

## Artifacts and persistence

By default, artifacts are written to the `artifacts/` directory:

- `artifacts/kmer_vectorizer.json`: saved vectorizer state (type plus parameters,
  e.g. k-mer vocabulary or transformer vectorizer configuration (DNABERT2/GenomeOcean))
- `artifacts/asv_vectors/{asv_id}.npy`: ASV-level vectors (one per unique ASV, reused across samples)
- `artifacts/sample_vectors/{sample_id}.npy`: sample-level vectors
- `artifacts/sample_index.faiss`: persisted FAISS index of sample vectors
- `artifacts/index_mapping.json`: ordered list mapping FAISS row indices to `sample_id`

 By default, the DuckDB database is created at the repo root as
 `reference.duckdb` with tables:

 - `samples(sample_id, metadata_id, sample_vector_path)`
 - `asvs(asv_id, sequence, abundance, sample_id, vector_path)`
 - `metadata(metadata_id, json_blob)`
 - `queries(query_id, timestamp, input_file, results_file)`

 When you pass a custom `db_path`, the same schema is created in the specified
 file.

## API usage

Start the API server:

```
uvicorn microbiome_tool.api.app:app --reload
```

Compare a sample by posting a TSV file:

```
curl -X POST \
  -F "file=@data/asv_table.tsv" \
  "http://127.0.0.1:8000/compare?top_n=5"
```

Optionally, you can pass a `db_path` query parameter to associate the request
with a particular DuckDB database file:

```
"http://127.0.0.1:8000/compare?top_n=5&db_path=my_reference.duckdb"
```

Response example:

```
{
  "results": [
    {"rank": 1, "sample_index": 0, "sample_id": "S1", "distance": 0.1234},
    {"rank": 2, "sample_index": 1, "sample_id": "S2", "distance": 0.2345}
  ]
}
```

## Programmatic usage (Python)

```python
from microbiome_tool.pipeline.ingestion import load_reference_data
from microbiome_tool.pipeline.comparison import compare_query_sample

# Build reference artifacts using the default k-mer vectorizer
res = load_reference_data(
    asv_table="data/asv_table.tsv",
    metadata="data/meta.tsv",
    k=5,
    artifacts_dir="artifacts",
)
print(res)  # {"samples": 3, "asvs": 9}

# Build reference artifacts from a wide ASV table (sequence x samples)
res_wide = load_reference_data(
    asv_table="data/asv_table_wide.tsv",
    metadata="data/meta.tsv",
    artifacts_dir="artifacts_wide",
    asv_table_format="wide",
)

# Compare a query against the FAISS index
results = compare_query_sample(
    query_table="data/asv_table.tsv",
    top_n=3,
    artifacts_dir="artifacts",
)
for r in results:
    print(r)

# Alternatively, build reference artifacts with DNABERT2 embeddings
res_dnabert2 = load_reference_data(
    asv_table="data/asv_table.tsv",
    metadata="data/meta.tsv",
    artifacts_dir="artifacts",
    vectorizer_type="dnabert2",
    dnabert2_model_name="zhihan1996/DNABERT-2-117M",
)

# Or build reference artifacts with GenomeOcean-100M embeddings
res_genomeocean = load_reference_data(
    asv_table="data/asv_table.tsv",
    metadata="data/meta.tsv",
    artifacts_dir="artifacts_genomeocean",
    vectorizer_type="genomeocean",
    genomeocean_model_name="DOEJGI/GenomeOcean-100M",
    genomeocean_model_max_length=1024,
)
```

To use a different DuckDB database file from Python, pass the optional
`db_path` argument when ingesting and (optionally) when comparing:

```python
res = load_reference_data(
    asv_table="data/asv_table.tsv",
    metadata="data/meta.tsv",
    artifacts_dir="artifacts_alt",
    db_path="my_reference.duckdb",
)

results = compare_query_sample(
    query_table="data/asv_table.tsv",
    top_n=3,
    artifacts_dir="artifacts_alt",
    db_path="my_reference.duckdb",
)
```

## Reproducibility and resetting

- By default, the `KmerVectorizer` builds and persists a sorted k-mer vocabulary so
  vector dimensions are consistent between ingestion and query. Alternative
  transformer-based vectorizers such as DNABERT2 or GenomeOcean-100M can also be
  used when ingesting.
- Vectors are L2-normalized by default.
- To reset a given reference, remove its artifacts directory (for example,
  `artifacts/` or your custom `artifacts_dir`) and the corresponding DuckDB
  file (for example, `reference.duckdb` or your custom `db_path`), then re-run
  ingestion.

## Extensibility

- Vectorizers are pluggable. To add a new vectorizer:
  - Subclass `Vectorizer` and implement `fit()` and `transform()`.
  - Provide JSON persistence helpers similar to the `KmerVectorizer`.
  - Update the ingestion/comparison steps to use your new vectorizer.

## Troubleshooting

- "No FAISS index found" or "No saved vectorizer found": run the ingestion step first.
- Dimension mismatch after changing `k` or vectorizer: delete `artifacts/` and re-run ingestion.
- Ensure your ASV table columns match the required schema.
