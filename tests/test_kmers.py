"""Unit tests for bac_probes.kmers."""

import pytest
from bac_probes.kmers import (
    reverse_complement,
    gc_content,
    _has_long_homopolymer,
    extract_kmers,
    score_conservation,
    score_offtarget_exact,
)


# ── reverse_complement ────────────────────────────────────────────────────────

def test_rc_basic():
    assert reverse_complement("ACGT") == "ACGT"

def test_rc_all_a():
    assert reverse_complement("AAAA") == "TTTT"

def test_rc_asymmetric():
    assert reverse_complement("ATCG") == "CGAT"

def test_rc_longer():
    seq = "GGTTTAATTCGATGATACGCGAGGAACCTTAC"
    rc  = "GTAAGGTTCCTCGCGTATCATCGAATTAAACC"
    assert reverse_complement(seq) == rc
    assert reverse_complement(rc) == seq


# ── gc_content ────────────────────────────────────────────────────────────────

def test_gc_all_at():
    assert gc_content("ATATAT") == 0.0

def test_gc_all_gc():
    assert gc_content("GCGCGC") == 1.0

def test_gc_half():
    assert gc_content("AACCGGTT") == pytest.approx(0.5)

def test_gc_32mer():
    # T-A-A-G-T-C-A-G-T-T-G-T-G-A-A-A-G-T-T-T-G-C-G-G-C-T-C-A-A-C-C-G
    # GC count = 15 (9 G + 6 C)
    kmer = "TAAGTCAGTTGTGAAAGTTTGCGGCTCAACCG"
    assert gc_content(kmer) == pytest.approx(15 / 32)


# ── _has_long_homopolymer ─────────────────────────────────────────────────────

def test_homopolymer_below_threshold():
    assert not _has_long_homopolymer("AAAACCCC", max_run=5)

def test_homopolymer_at_threshold():
    assert _has_long_homopolymer("AAAAACGT", max_run=5)

def test_homopolymer_single_base():
    assert not _has_long_homopolymer("ACGT", max_run=2)

def test_homopolymer_run_of_four():
    assert _has_long_homopolymer("ACGTTTTACGT", max_run=4)

def test_homopolymer_run_of_three_ok():
    assert not _has_long_homopolymer("ACGTTTACGT", max_run=4)


# ── extract_kmers ─────────────────────────────────────────────────────────────

def test_extract_kmers_basic():
    seq = "ACGTACGT"
    kmers = extract_kmers(seq, k=4, skip_ambiguous=True, gc_range=(0.0, 1.0), max_homopolymer=10)
    assert "ACGT" in kmers
    assert "CGTA" in kmers
    assert len(kmers) <= 5

def test_extract_kmers_ambiguous_filtered():
    seq = "ACGTNACGT"
    kmers = extract_kmers(seq, k=4)
    assert all("N" not in k for k in kmers)

def test_extract_kmers_gc_filter():
    seq = "AAAAAAAAACGTCGT"
    kmers = extract_kmers(seq, k=4, gc_range=(0.5, 1.0))
    for k in kmers:
        assert gc_content(k) >= 0.5

def test_extract_kmers_homopolymer_filter():
    seq = "AAAAACGTACGT"
    kmers = extract_kmers(seq, k=5, max_homopolymer=5)
    assert "AAAAA" not in kmers

def test_extract_kmers_returns_set():
    # Even if kmer appears multiple times, it's returned once
    seq = "ACGTACGTACGT"
    kmers = extract_kmers(seq, k=4)
    assert isinstance(kmers, set)

def test_extract_kmers_too_short():
    assert extract_kmers("ACG", k=4) == set()


# ── score_conservation ────────────────────────────────────────────────────────

def test_conservation_perfect():
    seqs = ["ACGTACGTACGT", "ACGTACGTACGT", "ACGTACGTACGT"]
    cons, n = score_conservation(seqs, k=4, min_conservation=1.0, gc_range=(0.0, 1.0))
    assert n == 3
    assert "ACGT" in cons
    assert cons["ACGT"] == pytest.approx(1.0)

def test_conservation_partial():
    seqs = ["ACGTTTTT", "ACGCCCCC", "ACGAAAAA"]
    cons, n = score_conservation(seqs, k=4, min_conservation=0.33, gc_range=(0.0, 1.0))
    assert n == 3
    # "ACGT" appears in 1/3, below 0.33 threshold
    # "ACGC" in 1/3, "ACGA" in 1/3 — all at the threshold
    assert cons.get("ACGT", 0) == pytest.approx(1/3, abs=1e-9)

def test_conservation_threshold_filter():
    seqs = ["ACGTTTTT", "ACGTCCCC", "AAAAGGGG"]
    # "ACGT" appears in 2/3 ≈ 0.667; below 0.70 threshold, above 0.60
    cons_60, _ = score_conservation(seqs, k=4, min_conservation=0.60, gc_range=(0.0, 1.0))
    cons_70, _ = score_conservation(seqs, k=4, min_conservation=0.70, gc_range=(0.0, 1.0))
    assert "ACGT" in cons_60
    assert "ACGT" not in cons_70

def test_conservation_empty_sequences():
    cons, n = score_conservation([], k=4)
    assert n == 0
    assert cons == {}

def test_conservation_respects_gc_filter():
    # "AAAA" has GC=0, should be filtered by gc_range
    seqs = ["AAAA" * 10, "AAAA" * 10]
    cons, n = score_conservation(seqs, k=4, min_conservation=0.5, gc_range=(0.35, 0.65))
    assert "AAAA" not in cons


# ── score_offtarget_exact ─────────────────────────────────────────────────────

def test_offtarget_exact_no_hits():
    candidates = {"ACGTACGTACGTACGTACGTACGTACGTACGT"}
    bg = ["TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"]
    counts, n = score_offtarget_exact(candidates, bg, k=32)
    assert n == 2
    assert counts["ACGTACGTACGTACGTACGTACGTACGTACGT"] == 0

def test_offtarget_exact_one_hit():
    kmer = "ACGT" * 8  # 32-mer
    candidates = {kmer}
    bg = [kmer, "T" * 32]
    counts, n = score_offtarget_exact(candidates, bg, k=32)
    assert counts[kmer] == 1
    assert n == 2

def test_offtarget_exact_multiple_candidates():
    k1 = "A" * 16 + "C" * 16
    k2 = "G" * 16 + "T" * 16
    candidates = {k1, k2}
    bg = [k1]  # only k1 is in background
    counts, n = score_offtarget_exact(candidates, bg, k=32)
    assert counts[k1] == 1
    assert counts[k2] == 0
