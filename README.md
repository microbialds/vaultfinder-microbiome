# Microbiome Tool

A research prototype for comparing microbiome samples based on amplicon sequencing (ASV tables).

The tool ingests ASV tables (DADA2 format), stores data in DuckDB, generates vectors (ASV + sample level),
indexes sample vectors with FAISS, and allows sample-level comparison via CLI or a minimal FastAPI service.

## Requirements

- Python >= 3.10
- See `requirements.txt` for dependencies (FAISS CPU, DuckDB, FastAPI, etc.)

## Installation

```bash
pip install -e .
```

## Usage

### Data formats

- ASV table (TSV), either:
  - Long format (default): columns `sample_id`, `asv_id`, `sequence`, `abundance`
  - Wide format: first column is the DNA sequence, subsequent columns are sample IDs with abundances
- Optional metadata (TSV) with at least `sample_id` and any additional fields

### Quickstart (CLI)

1) Build reference artifacts (vectorizer, vectors, FAISS index, mapping):

```bash
load-reference --asv_table data/asv_table.tsv --metadata data/meta.tsv
```

To ingest a wide-format ASV table (sequence x samples):

```bash
load-reference --asv_table data/asv_table_wide.tsv --metadata data/meta.tsv --asv_table_format wide
```

2) Compare a query sample against the FAISS index:

```bash
compare-samples --query data/asv_table.tsv --top_n 5
```

By default, both commands use `reference.duckdb` in the repo root as the DuckDB
database. To write to a different database file when ingesting, add
`--db_path my_reference.duckdb` to `load-reference`. When comparing, you can
also pass `--artifacts_dir` and `--db_path` if you manage multiple, separate
references.

To use DNABERT2 embeddings instead of k-mers when building the reference, pass the
`dnabert2` vectorizer type and (optionally) a model name:

```bash
load-reference \
  --asv_table data/asv_table.tsv \
  --metadata data/meta.tsv \
  --vectorizer_type dnabert2 \
  --dnabert2_model_name zhihan1996/DNABERT-2-117M
```

To use GenomeOcean-100M (Mistral-based) embeddings instead of k-mers, pass the
`genomeocean` vectorizer type and (optionally) a model name and maximum token
length. The default model is `DOEJGI/GenomeOcean-100M`:

```bash
load-reference \
  --asv_table data/asv_table.tsv \
  --metadata data/meta.tsv \
  --vectorizer_type genomeocean \
  --genomeocean_model_name DOEJGI/GenomeOcean-100M \
  --genomeocean_model_max_length 1024 \
  --device cuda
```

### API

Start the API server:

```bash
uvicorn microbiome_tool.api.app:app --reload
```

Compare via HTTP by posting a TSV file:

```bash
curl -X POST \
  -F "file=@data/asv_table.tsv" \
  "http://127.0.0.1:8000/compare?top_n=5"
```

### Programmatic (Python)

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

To target a different DuckDB database file from Python, pass the optional
`db_path` argument when ingesting and (optionally) when comparing:

```python
load_reference_data(
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

### Artifacts & persistence

By default, artifacts are written to `artifacts/`:

- `kmer_vectorizer.json`: saved vectorizer state (type plus parameters, e.g. k-mer
  vocabulary or transformer vectorizer configuration (DNABERT2/GenomeOcean))
- `asv_vectors/{asv_id}.npy`: ASV-level vectors (one per unique ASV, reused across samples)
- `sample_vectors/{sample_id}.npy`: sample-level vectors
- `sample_index.faiss`: FAISS index of sample vectors
- `index_mapping.json`: mapping from FAISS row indices to `sample_id`

By default, the DuckDB database is `reference.duckdb` with tables
`samples`, `asvs`, `metadata`, `queries`. When using a custom `db_path`, the
same schema is created in the specified file.

To reset a given reference, delete its artifacts directory (for example,
`artifacts/` or your custom `artifacts_dir`) and the corresponding DuckDB file
(for example, `reference.duckdb` or your custom `db_path`), then re-run
ingestion.

## Layout

- `microbiome_tool/`: core package (DuckDB, FAISS, vectorizers, pipeline, API)
- `cli/`: CLI entry points
- `data/`: example inputs
- `tests/`: pytest scaffold
- `docs/`: additional documentation (`USAGE.md`, `ARCHITECTURE.md`)

## Documentation

- Usage guide: `docs/USAGE.md`
- Architecture: `docs/ARCHITECTURE.md`

## License

MIT License
