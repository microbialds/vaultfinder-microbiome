import hashlib
import pandas as pd

def read_asv_table(path: str) -> pd.DataFrame:
    """Read a DADA2-like ASV table (TSV) with columns: sample_id, asv_id, sequence, abundance."""
    return pd.read_csv(path, sep="\t")


def read_asv_table_wide(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return df

    # Ensure the first column is treated as the sequence column; leave the
    # remaining columns as-is (one column per sample with abundances).
    seq_col = df.columns[0]
    df[seq_col] = df[seq_col].astype(str)
    return df


def hash_sequence(sequence: str) -> str:
    return hashlib.md5(sequence.encode("utf-8")).hexdigest()
