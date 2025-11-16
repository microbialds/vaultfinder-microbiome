import numpy as np
from collections import Counter
from typing import List, Optional
import json
import os

class Vectorizer:
    def fit(self, sequences):
        """Optional fit for vectorizers that need training."""
        return self

    def transform(self, sequences, abundances):
        """Return a vector representation for the (sequences, abundances)."""
        raise NotImplementedError

class KmerVectorizer(Vectorizer):
    def __init__(self, k=5, vocab: Optional[List[str]] = None):
        self.k = k
        # When set via fit() or provided explicitly, ensures fixed-dimension outputs
        self.vocab_: Optional[List[str]] = list(vocab) if vocab is not None else None

    def fit(self, sequences: List[str]):
        """Build a deterministic k-mer vocabulary from sequences.

        The vocabulary is the sorted list of all unique k-mers observed across the
        provided sequences. Using a fixed vocabulary guarantees that all vectors have
        the same dimensionality and consistent ordering across runs.
        """
        vocab_set = set()
        for seq in sequences:
            if not isinstance(seq, str):
                continue
            if len(seq) < self.k:
                continue
            for i in range(len(seq) - self.k + 1):
                vocab_set.add(seq[i:i + self.k])
        self.vocab_ = sorted(vocab_set)
        return self

    def _kmerize(self, seq):
        return [seq[i:i+self.k] for i in range(len(seq) - self.k + 1)]

    def transform(self, sequences, abundances, normalize: bool = True):
        counts = Counter()
        for seq, ab in zip(sequences, abundances):
            if not isinstance(seq, str) or ab is None:
                continue
            kmers = self._kmerize(seq)
            for kmer in kmers:
                counts[kmer] += ab

        # If a vocabulary is set, produce a fixed-length vector in that order
        if getattr(self, "vocab_", None):
            vec = np.zeros(len(self.vocab_), dtype=float)
            if counts:
                # Map counts into the fixed vocab order
                vocab_index = {k: i for i, k in enumerate(self.vocab_)}
                for kmer, c in counts.items():
                    idx = vocab_index.get(kmer)
                    if idx is not None:
                        vec[idx] = c
            if normalize and np.linalg.norm(vec) > 0:
                vec = vec / np.linalg.norm(vec)
            return vec

        # Backward-compatible behavior when no vocabulary is defined: variable-length
        keys = sorted(counts.keys())
        vec = np.array([counts[k] for k in keys], dtype=float)
        if normalize and np.linalg.norm(vec) > 0:
            vec = vec / np.linalg.norm(vec)
        return vec

def save_kmer_vectorizer(vec: KmerVectorizer, path: str):
    """Persist KmerVectorizer parameters (k, vocab_) as JSON."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        "type": "kmer",
        "k": vec.k,
        "vocab": getattr(vec, "vocab_", None) or [],
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f)

def load_kmer_vectorizer(path: str) -> KmerVectorizer:
    """Load a KmerVectorizer from a JSON file produced by save_kmer_vectorizer."""
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if payload.get("type") != "kmer":
        raise ValueError("Unsupported vectorizer type in saved file")
    k = int(payload.get("k", 5))
    vocab = payload.get("vocab", [])
    return KmerVectorizer(k=k, vocab=vocab)


class DNABERT2Vectorizer(Vectorizer):
    def __init__(
        self,
        model_name: str = "zhihan1996/DNABERT-2-117M",
        device: Optional[str] = None,
        batch_size: int = 8,
        pooling: str = "cls",
    ):
        self.model_name = model_name
        self.device = device
        self.batch_size = batch_size
        self.pooling = pooling
        self._tokenizer = None
        self._model = None

    def _ensure_model(self):
        if self._model is not None and self._tokenizer is not None:
            return
        try:
            from transformers import AutoTokenizer, AutoModel
            import torch  # noqa: F401  # used for device placement
        except ImportError as e:
            raise ImportError(
                "DNABERT2Vectorizer requires the 'transformers' and 'torch' packages "
                "to be installed."
            ) from e

        self._tokenizer = AutoTokenizer.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self._model = AutoModel.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self._model.eval()

        if self.device is None:
            if torch.cuda.is_available():
                self.device = "cuda"
            else:
                self.device = "cpu"

        if self.device is not None:
            self._model.to(self.device)

    def fit(self, sequences):
        return self

    def transform(self, sequences, abundances, normalize: bool = True):
        self._ensure_model()

        import numpy as _np
        import torch

        filtered_seqs = []
        weights = []
        for seq, ab in zip(sequences, abundances):
            if not isinstance(seq, str) or ab is None:
                continue
            filtered_seqs.append(seq)
            weights.append(float(ab))

        if not filtered_seqs:
            hidden_size = getattr(getattr(self._model, "config", None), "hidden_size", 0)
            if hidden_size <= 0:
                return np.zeros(0, dtype=float)
            return np.zeros(hidden_size, dtype=float)

        batch_size = max(1, int(self.batch_size))
        all_embs = []
        all_weights = []

        for i in range(0, len(filtered_seqs), batch_size):
            batch_seqs = filtered_seqs[i : i + batch_size]
            batch_weights = weights[i : i + batch_size]

            inputs = self._tokenizer(
                batch_seqs,
                return_tensors="pt",
                padding=True,
                truncation=True,
            )
            if self.device is not None:
                inputs = {k: v.to(self.device) for k, v in inputs.items()}

            with torch.no_grad():
                outputs = self._model(**inputs)

            if isinstance(outputs, tuple):
                last_hidden = outputs[0]
            else:
                last_hidden = outputs.last_hidden_state

            if self.pooling == "mean":
                mask = inputs["attention_mask"].unsqueeze(-1)
                summed = (last_hidden * mask).sum(dim=1)
                counts = mask.sum(dim=1).clamp(min=1)
                batch_emb = summed / counts
            else:
                batch_emb = last_hidden[:, 0, :]

            batch_emb = batch_emb.cpu().numpy()
            all_embs.append(batch_emb)
            all_weights.extend(batch_weights)

        embs = _np.concatenate(all_embs, axis=0)
        w = _np.asarray(all_weights, dtype=float)

        total = w.sum()
        if total > 0:
            vec = (embs * w[:, None]).sum(axis=0) / total
        else:
            vec = embs.mean(axis=0)

        if normalize and _np.linalg.norm(vec) > 0:
            vec = vec / _np.linalg.norm(vec)

        return vec.astype(float)


def save_dnabert2_vectorizer(vec: DNABERT2Vectorizer, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        "type": "dnabert2",
        "model_name": vec.model_name,
        "pooling": vec.pooling,
        "batch_size": vec.batch_size,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f)


def load_dnabert2_vectorizer(path: str) -> DNABERT2Vectorizer:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if payload.get("type") != "dnabert2":
        raise ValueError("Unsupported vectorizer type in saved file")
    model_name = payload.get("model_name", "zhihan1996/DNABERT-2-117M")
    pooling = payload.get("pooling", "cls")
    batch_size = int(payload.get("batch_size", 8))
    return DNABERT2Vectorizer(
        model_name=model_name,
        pooling=pooling,
        batch_size=batch_size,
    )


class GenomeOceanVectorizer(Vectorizer):
    def __init__(
        self,
        model_name: str = "DOEJGI/GenomeOcean-100M",
        device: Optional[str] = None,
        batch_size: int = 8,
        pooling: str = "mean",
        model_max_length: Optional[int] = None,
    ):
        self.model_name = model_name
        self.device = device
        self.batch_size = batch_size
        self.pooling = pooling
        self.model_max_length = model_max_length
        self._tokenizer = None
        self._model = None

    def _ensure_model(self):
        if self._model is not None and self._tokenizer is not None:
            return
        try:
            from transformers import AutoTokenizer, AutoModelForCausalLM
            import torch  # noqa: F401  # used for device placement
        except ImportError as e:
            raise ImportError(
                "GenomeOceanVectorizer requires the 'transformers' and 'torch' packages "
                "to be installed."
            ) from e

        self._tokenizer = AutoTokenizer.from_pretrained(
            self.model_name,
            trust_remote_code=True,
            padding_side="left",
        )
        self._model = AutoModelForCausalLM.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self._model.eval()

        # Disable KV caching to avoid incompatibilities with DynamicCache/
        # get_usable_length across different transformers versions.
        cfg = getattr(self._model, "config", None)
        if cfg is not None and hasattr(cfg, "use_cache"):
            cfg.use_cache = False

        if self.device is None:
            if torch.cuda.is_available():
                self.device = "cuda"
            else:
                self.device = "cpu"

        if self.device is not None:
            self._model.to(self.device)

    def fit(self, sequences):
        return self

    def transform(self, sequences, abundances, normalize: bool = True):
        self._ensure_model()

        import numpy as _np
        import torch

        filtered_seqs = []
        weights = []
        for seq, ab in zip(sequences, abundances):
            if not isinstance(seq, str) or ab is None:
                continue
            filtered_seqs.append(seq)
            weights.append(float(ab))

        if not filtered_seqs:
            hidden_size = getattr(getattr(self._model, "config", None), "hidden_size", 0)
            if hidden_size <= 0:
                return np.zeros(0, dtype=float)
            return np.zeros(hidden_size, dtype=float)

        batch_size = max(1, int(self.batch_size))
        all_embs = []
        all_weights = []

        for i in range(0, len(filtered_seqs), batch_size):
            batch_seqs = filtered_seqs[i : i + batch_size]
            batch_weights = weights[i : i + batch_size]

            inputs = self._tokenizer(
                batch_seqs,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=self.model_max_length,
            )
            # Some causal LM implementations (e.g., Mistral) do not accept
            # token_type_ids; remove them if present in the tokenizer output.
            if "token_type_ids" in inputs:
                inputs.pop("token_type_ids")
            if self.device is not None:
                inputs = {k: v.to(self.device) for k, v in inputs.items()}

            with torch.no_grad():
                outputs = self._model(
                    **inputs,
                    use_cache=False,
                    output_hidden_states=True,
                    return_dict=True,
                )

            hidden_states = outputs.hidden_states[-1]

            if self.pooling == "mean":
                mask = inputs["attention_mask"].unsqueeze(-1)
                summed = (hidden_states * mask).sum(dim=1)
                counts = mask.sum(dim=1).clamp(min=1)
                batch_emb = summed / counts
            else:
                # Last-token pooling, respecting left padding
                mask = inputs["attention_mask"]
                lengths = mask.sum(dim=1) - 1
                lengths = lengths.clamp(min=0)
                batch_indices = torch.arange(
                    hidden_states.size(0), device=hidden_states.device
                )
                batch_emb = hidden_states[batch_indices, lengths, :]

            batch_emb = batch_emb.cpu().numpy()
            all_embs.append(batch_emb)
            all_weights.extend(batch_weights)

        embs = _np.concatenate(all_embs, axis=0)
        w = _np.asarray(all_weights, dtype=float)

        total = w.sum()
        if total > 0:
            vec = (embs * w[:, None]).sum(axis=0) / total
        else:
            vec = embs.mean(axis=0)

        if normalize and _np.linalg.norm(vec) > 0:
            vec = vec / _np.linalg.norm(vec)

        return vec.astype(float)


def save_genomeocean_vectorizer(vec: GenomeOceanVectorizer, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload = {
        "type": "genomeocean",
        "model_name": vec.model_name,
        "pooling": vec.pooling,
        "batch_size": vec.batch_size,
        "model_max_length": vec.model_max_length,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f)


def load_genomeocean_vectorizer(path: str) -> GenomeOceanVectorizer:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if payload.get("type") != "genomeocean":
        raise ValueError("Unsupported vectorizer type in saved file")
    model_name = payload.get("model_name", "DOEJGI/GenomeOcean-100M")
    pooling = payload.get("pooling", "mean")
    batch_size = int(payload.get("batch_size", 8))
    model_max_length = payload.get("model_max_length", None)
    return GenomeOceanVectorizer(
        model_name=model_name,
        pooling=pooling,
        batch_size=batch_size,
        model_max_length=model_max_length,
    )


def save_vectorizer(vec: Vectorizer, path: str):
    if isinstance(vec, KmerVectorizer):
        save_kmer_vectorizer(vec, path)
        return
    if isinstance(vec, DNABERT2Vectorizer):
        save_dnabert2_vectorizer(vec, path)
        return
    if isinstance(vec, GenomeOceanVectorizer):
        save_genomeocean_vectorizer(vec, path)
        return
    raise ValueError("Unsupported vectorizer instance for saving")


def load_vectorizer(path: str) -> Vectorizer:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    vtype = payload.get("type")
    if vtype == "kmer":
        return load_kmer_vectorizer(path)
    if vtype == "dnabert2":
        return load_dnabert2_vectorizer(path)
    if vtype == "genomeocean":
        return load_genomeocean_vectorizer(path)
    raise ValueError("Unsupported vectorizer type in saved file")
