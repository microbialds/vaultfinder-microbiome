import click
from microbiome_tool.pipeline.comparison import compare_query_sample

@click.command()
@click.option("--query", required=True, help="ASV table for query sample")
@click.option("--top_n", default=5, help="Number of top matches")
@click.option(
    "--artifacts_dir",
    default="artifacts",
    show_default=True,
    help="Directory where artifacts (vectorizer, FAISS, mapping) are stored",
)
@click.option(
    "--db_path",
    required=False,
    help="Optional DuckDB database path (defaults to reference.duckdb)",
)
def compare_samples(query, top_n, artifacts_dir, db_path):
    results = compare_query_sample(
        query_table=query,
        top_n=top_n,
        artifacts_dir=artifacts_dir,
        db_path=db_path,
    )
    print("Top matches:")
    for r in results:
        sid = r.get('sample_id') or 'N/A'
        print(f"{r['rank']}. SampleIdx {r['sample_index']} (ID: {sid}) | Distance: {r['distance']:.4f}")

if __name__ == "__main__":
    compare_samples()
