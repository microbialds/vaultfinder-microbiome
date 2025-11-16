import click
from microbiome_tool.pipeline.ingestion import load_reference_data

@click.command()
@click.option("--asv_table", required=True, help="ASV table in DADA2 format")
@click.option(
    "--asv_table_format",
    type=click.Choice(["long", "wide"], case_sensitive=False),
    default="long",
    show_default=True,
    help="ASV table format: 'long' (sample_id/asv_id/sequence/abundance) or 'wide' (sequence x samples)",
)
@click.option("--metadata", required=False, help="Optional metadata file")
@click.option("--k", type=int, default=5, show_default=True, help="k-mer size for KmerVectorizer")
@click.option(
    "--artifacts_dir",
    default="artifacts",
    show_default=True,
    help="Directory where artifacts (vectors, FAISS, mapping) will be written",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Suppress progress output",
)
@click.option(
    "--vectorizer_type",
    type=click.Choice(["kmer", "dnabert2", "genomeocean"], case_sensitive=False),
    default="kmer",
    show_default=True,
    help="Vectorizer type to use when building sample vectors",
)
@click.option(
    "--dnabert2_model_name",
    required=False,
    help="Optional DNABERT2 Hugging Face model name",
)
@click.option(
    "--genomeocean_model_name",
    required=False,
    help="Optional GenomeOcean Hugging Face model name (defaults to DOEJGI/GenomeOcean-100M)",
)
@click.option(
    "--genomeocean_model_max_length",
    type=int,
    required=False,
    help="Optional maximum token length for GenomeOcean tokenizer (passed as max_length)",
)
@click.option(
    "--device",
    type=click.Choice(["auto", "cpu", "cuda"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="Device to use for transformer-based vectorizers (DNABERT2 or GenomeOcean)",
)
@click.option(
    "--db_path",
    required=False,
    help="Optional DuckDB database path (defaults to reference.duckdb)",
)
def load_reference(
    asv_table,
    asv_table_format,
    metadata,
    k,
    artifacts_dir,
    quiet,
    vectorizer_type,
    dnabert2_model_name,
    genomeocean_model_name,
    genomeocean_model_max_length,
    device,
    db_path,
):
    # Default to verbose unless quiet is explicitly requested
    verbose = not quiet
    result = load_reference_data(
        asv_table,
        metadata,
        k=k,
        artifacts_dir=artifacts_dir,
        verbose=verbose,
        vectorizer_type=vectorizer_type,
        dnabert2_model_name=dnabert2_model_name,
        dnabert2_device=device,
        genomeocean_model_name=genomeocean_model_name,
        genomeocean_model_max_length=genomeocean_model_max_length,
        db_path=db_path,
        asv_table_format=asv_table_format,
    )
    if verbose and isinstance(result, dict):
        click.echo(
            f"Completed ingestion: samples={result.get('samples', 0)}, asvs={result.get('asvs', 0)}"
        )

if __name__ == "__main__":
    load_reference()
