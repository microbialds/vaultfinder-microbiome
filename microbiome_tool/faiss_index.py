import faiss
import numpy as np
import os

INDEX_FILE = "sample_index.faiss"

def create_index(dim):
    """Create a simple L2 FAISS index with given dimensionality."""
    index = faiss.IndexFlatL2(dim)
    return index

def save_index(index, path=INDEX_FILE):
    faiss.write_index(index, path)

def load_index(path=INDEX_FILE):
    if os.path.exists(path):
        return faiss.read_index(path)
    return None

def add_to_index(index, vectors):
    index.add(vectors)  # vectors: (n, dim) float32
    return index

def search_index(index, query_vec, top_n=5):
    distances, indices = index.search(query_vec, top_n)
    return distances, indices
