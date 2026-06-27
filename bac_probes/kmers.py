"""K-mer extraction, conservation scoring, and exact off-target counting."""

from collections import Counter
from typing import Iterable


VALID_BASES = frozenset("ACGT")


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def gc_content(kmer: str) -> float:
    return sum(c in "GC" for c in kmer) / len(kmer)


def _has_long_homopolymer(kmer: str, max_run: int = 5) -> bool:
    """Return True if kmer contains a run of the same base >= max_run."""
    run_base, run_len = "", 0
    for base in kmer:
        if base == run_base:
            run_len += 1
            if run_len >= max_run:
                return True
        else:
            run_base, run_len = base, 1
    return False


def extract_kmers(
    sequence: str,
    k: int,
    skip_ambiguous: bool = True,
    gc_range: tuple[float, float] = (0.0, 1.0),
    max_homopolymer: int = 5,
) -> set[str]:
    """
    Return the set of unique k-mers in a sequence that pass quality filters.

    Filters applied (when enabled):
    - skip_ambiguous: drop k-mers containing any base not in ACGT
    - gc_range: drop k-mers whose GC fraction is outside [min, max]
    - max_homopolymer: drop k-mers with homopolymer runs >= this length
    """
    kmers: set[str] = set()
    n = len(sequence)
    for i in range(n - k + 1):
        kmer = sequence[i : i + k]
        if skip_ambiguous and not VALID_BASES.issuperset(kmer):
            continue
        gc = gc_content(kmer)
        if not (gc_range[0] <= gc <= gc_range[1]):
            continue
        if _has_long_homopolymer(kmer, max_homopolymer):
            continue
        kmers.add(kmer)
    return kmers


def score_conservation(
    sequences: Iterable[str],
    k: int = 32,
    min_conservation: float = 0.80,
    gc_range: tuple[float, float] = (0.35, 0.65),
    max_homopolymer: int = 5,
) -> tuple[dict[str, float], int]:
    """
    Scan target sequences for k-mers and return per-kmer conservation fractions.

    A k-mer's conservation is the fraction of target sequences that contain it
    (at least once, anywhere in the sequence).

    Returns:
        conservation: {kmer: fraction} for k-mers with fraction >= min_conservation
        n_seqs: number of sequences processed
    """
    kmer_counts: Counter = Counter()
    n_seqs = 0

    for seq in sequences:
        kmers = extract_kmers(seq, k, gc_range=gc_range, max_homopolymer=max_homopolymer)
        kmer_counts.update(kmers)
        n_seqs += 1

    if n_seqs == 0:
        return {}, 0

    conservation = {
        kmer: count / n_seqs
        for kmer, count in kmer_counts.items()
        if count / n_seqs >= min_conservation
    }
    return conservation, n_seqs


def score_offtarget_exact(
    candidate_kmers: set[str],
    background_sequences: Iterable[str],
    k: int = 32,
) -> tuple[dict[str, int], int]:
    """
    Count how many background (non-target) sequences contain each candidate k-mer.

    Uses streaming to avoid storing background sequences in memory.

    Returns:
        offtarget_counts: {kmer: count_of_bg_seqs_containing_it}
        n_bg_seqs: total number of background sequences scanned
    """
    offtarget_counts: dict[str, int] = {kmer: 0 for kmer in candidate_kmers}
    n_bg_seqs = 0

    for seq in background_sequences:
        # Extract k-mers with no extra filters so we catch everything
        seq_kmers: set[str] = set()
        n = len(seq)
        for i in range(n - k + 1):
            kmer = seq[i : i + k]
            if VALID_BASES.issuperset(kmer):
                seq_kmers.add(kmer)

        for kmer in candidate_kmers & seq_kmers:
            offtarget_counts[kmer] += 1

        n_bg_seqs += 1

    return offtarget_counts, n_bg_seqs
