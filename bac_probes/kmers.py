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
    Handles variable-length candidates (from merge_overlapping_kmers) via substring search.

    Returns:
        offtarget_counts: {kmer: count_of_bg_seqs_containing_it}
        n_bg_seqs: total number of background sequences scanned
    """
    offtarget_counts: dict[str, int] = {kmer: 0 for kmer in candidate_kmers}
    n_bg_seqs = 0
    variable_length = any(len(c) != k for c in candidate_kmers)

    for seq in background_sequences:
        if variable_length:
            for kmer in candidate_kmers:
                if kmer in seq:
                    offtarget_counts[kmer] += 1
        else:
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


def score_conservation_with_coverage(
    sequences: list[str],
    k: int = 32,
    min_conservation: float = 0.20,
    gc_range: tuple[float, float] = (0.35, 0.65),
    max_homopolymer: int = 5,
) -> tuple[dict[str, set[int]], int]:
    """
    Like score_conservation but records which sequence indices contain each k-mer.

    Returns:
        coverage_sets: {kmer: set_of_seq_indices}
        n_seqs: number of sequences processed
    """
    kmer_to_seqs: dict[str, set[int]] = {}
    n_seqs = len(sequences)

    for i, seq in enumerate(sequences):
        for kmer in extract_kmers(seq, k, gc_range=gc_range, max_homopolymer=max_homopolymer):
            kmer_to_seqs.setdefault(kmer, set()).add(i)

    threshold = min_conservation * n_seqs
    return {k: s for k, s in kmer_to_seqs.items() if len(s) >= threshold}, n_seqs


def greedy_cocktail(
    coverage_sets: dict[str, set[int]],
    specificity: dict[str, float],
    n_total: int,
    n_probes: int = 8,
    coverage_target: float = 0.90,
    min_specificity: float = 0.0,
) -> list[str]:
    """
    Greedy set-cover: select up to n_probes k-mers that collectively cover
    >= coverage_target fraction of target sequences, ranked by
    (newly covered seqs) × specificity at each step.

    Returns ordered list of selected k-mers (best first).
    """
    candidates = {k: v for k, v in coverage_sets.items()
                  if specificity.get(k, 0.0) >= min_specificity}
    uncovered = set(range(n_total))
    selected: list[str] = []

    while len(selected) < n_probes and uncovered and candidates:
        best = max(
            candidates,
            key=lambda k: len(candidates[k] & uncovered) * specificity.get(k, 0.0),
        )
        newly = candidates[best] & uncovered
        if not newly:
            break
        selected.append(best)
        uncovered -= newly
        del candidates[best]
        if (n_total - len(uncovered)) / n_total >= coverage_target:
            break

    return selected


def merge_overlapping_kmers(kmer_set: set[str], k: int) -> list[str]:
    """
    Merge k-mers that overlap by exactly k-1 bases into longer contiguous sequences.

    Two k-mers A and B overlap if A[1:] == B[:k-1] (they are consecutive positions
    in the same 16S sequence). Chains of overlapping k-mers are assembled into a
    single longer sequence. Branches (ambiguous next k-mer) break the chain.
    """
    if not kmer_set:
        return []

    head_to_kmers: dict[str, list[str]] = {}
    for kmer in kmer_set:
        head_to_kmers.setdefault(kmer[:k - 1], []).append(kmer)

    all_tails = {kmer[1:] for kmer in kmer_set}
    start_kmers = [kmer for kmer in kmer_set if kmer[:k - 1] not in all_tails]

    merged: list[str] = []
    visited: set[str] = set()

    for start in start_kmers:
        if start in visited:
            continue
        chain = [start]
        visited.add(start)
        current = start
        while True:
            tail = current[1:]
            nexts = [n for n in head_to_kmers.get(tail, []) if n not in visited]
            if len(nexts) != 1:
                break
            nxt = nexts[0]
            chain.append(nxt)
            visited.add(nxt)
            current = nxt

        merged_seq = chain[0] + "".join(c[-1] for c in chain[1:])
        merged.append(merged_seq)

    for kmer in kmer_set:
        if kmer not in visited:
            merged.append(kmer)

    return merged


def rescore_merged_conservation(
    sequences: list[str],
    candidates: list[str],
    min_conservation: float = 0.80,
    gc_range: tuple[float, float] = (0.35, 0.65),
    max_homopolymer: int = 5,
) -> dict[str, float]:
    """
    Score conservation of variable-length sequences as fraction of target seqs containing them.
    Re-applies GC and homopolymer filters to the merged sequences.
    """
    n_seqs = len(sequences)
    if n_seqs == 0:
        return {}

    valid = [
        c for c in candidates
        if (gc_range[0] <= gc_content(c) <= gc_range[1])
        and not _has_long_homopolymer(c, max_homopolymer)
        and VALID_BASES.issuperset(c)
    ]

    counts = {c: 0 for c in valid}
    for seq in sequences:
        for c in valid:
            if c in seq:
                counts[c] += 1

    return {
        c: count / n_seqs
        for c, count in counts.items()
        if count / n_seqs >= min_conservation
    }
