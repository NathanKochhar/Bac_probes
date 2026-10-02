"""Unit tests for bac_probes.specificity."""

import pytest
from bac_probes.specificity import (
    mismatch_weight,
    _parse_btop,
    _hit_mismatch_positions,
    effective_binding_score,
    blast_batch_size,
    parse_blast_results,
)


# ── mismatch_weight ───────────────────────────────────────────────────────────

def test_weight_center_is_max():
    k = 32
    # Centre position (15 or 16) should be near the maximum
    w_center = mismatch_weight(15, k)
    w_end = mismatch_weight(0, k)
    assert w_center > w_end

def test_weight_ends_near_zero():
    assert mismatch_weight(0, 32) == pytest.approx(0.00239, abs=1e-4)
    assert mismatch_weight(31, 32) == pytest.approx(0.00239, abs=1e-4)

def test_weight_center_near_one():
    # sin(pi * 15.5 / 32)^2 ≈ sin(pi/2)^2 = 1 for the midpoint
    assert mismatch_weight(15, 32) == pytest.approx(1.0, abs=0.01)

def test_weight_symmetric():
    k = 32
    for i in range(k // 2):
        assert mismatch_weight(i, k) == pytest.approx(mismatch_weight(k - 1 - i, k), abs=1e-12)

def test_weight_quarter_position():
    # sin(π * 7.5 / 32)^2 ≈ 0.451; exact half-max falls between positions 7 and 8
    assert mismatch_weight(7, 32) == pytest.approx(0.451, abs=0.005)


# ── _parse_btop ───────────────────────────────────────────────────────────────

def test_btop_all_matches():
    assert _parse_btop("32") == []

def test_btop_single_mismatch():
    # 15 matches, mismatch at pos 15, 16 matches
    result = _parse_btop("15GA16")
    assert result == [15]

def test_btop_two_mismatches():
    result = _parse_btop("5AT10GA10")
    assert result == [5, 16]

def test_btop_gap_in_query():
    # "-G" = gap in query: subject advances, query does not. The bulge is
    # still penalised, at the query base following the insertion.
    result = _parse_btop("10-G10")
    assert result == [10]

def test_btop_qstart_offset():
    # Alignment starting at query base 3 (1-based) shifts positions by 2
    assert _parse_btop("4TC11AT12", qstart=3) == [6, 18]


# ── _hit_mismatch_positions ───────────────────────────────────────────────────

def test_hit_positions_full_alignment():
    assert _hit_mismatch_positions("15GA16", qstart=1, qend=32, qlen=32) == [15]

def test_hit_positions_include_unaligned_overhangs():
    # Probe bases 0-1 and 31 were not aligned; mismatches offset by qstart
    assert _hit_mismatch_positions("4TC11AT12", qstart=3, qend=31, qlen=32) == [0, 1, 31, 6, 18]

def test_btop_gap_in_subject():
    # "A-" = gap in subject: query advances by 1
    result = _parse_btop("10A-10")
    assert result == [10]

def test_btop_leading_mismatch():
    result = _parse_btop("AG30")
    assert result == [0]

def test_btop_trailing_mismatch():
    result = _parse_btop("31CT")
    assert result == [31]

def test_btop_no_match_prefix():
    # pure mismatch string
    result = _parse_btop("AGCT")
    assert result == [0, 1]


# ── effective_binding_score ───────────────────────────────────────────────────

def test_binding_perfect_match():
    # No mismatches → no protection → 0.0 (off-target will bind)
    assert effective_binding_score([], k=32) == 0.0

def test_binding_center_mismatch_maximal_protection():
    # Single central mismatch (pos 15) has weight ≈ 1.0 → score ≈ 1.0 (well protected)
    score = effective_binding_score([15], k=32)
    assert score == pytest.approx(1.0, abs=0.02)

def test_binding_end_mismatch_minimal_protection():
    # End mismatch (pos 0) has weight ≈ 0.002 → score ≈ 0.002 (minimal protection)
    score = effective_binding_score([0], k=32)
    assert score == pytest.approx(0.002, abs=0.01)

def test_binding_clamped_at_one():
    # Multiple heavy mismatches should not exceed 1.0
    score = effective_binding_score([14, 15, 16], k=32)
    assert score <= 1.0

def test_binding_two_end_mismatches():
    # Two terminal mismatches provide minimal combined protection (~0.005)
    score = effective_binding_score([0, 31], k=32)
    assert score < 0.01


# ── parse_blast_results ───────────────────────────────────────────────────────

K = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bp test kmer

def _make_blast_line(qseqid, stitle, pident, btop, length=32, qlen=32):
    return f"{qseqid}\t{stitle}\t{pident}\t{length}\t{qlen}\t{btop}"


def test_parse_on_target_exact_match():
    # stitle is the bare taxonomy string (as produced by -parse_seqids makeblastdb)
    stitle = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides fragilis"
    line = _make_blast_line("kmer_0", stitle, 100.0, "32")
    results = parse_blast_results(line, [K], "Bacteroides", "genus")

    r = results[K]
    assert r["blast_target_hits"] == 1
    assert r["blast_offtarget_hits"] == 0
    assert r["blast_specificity"] == 1.0
    assert r["blast_weighted_specificity"] == 1.0


def test_parse_off_target_exact_match():
    stitle = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;uncultured bacterium"
    line = _make_blast_line("kmer_0", stitle, 100.0, "32")
    results = parse_blast_results(line, [K], "Bacteroides", "genus")

    r = results[K]
    assert r["blast_target_hits"] == 0
    assert r["blast_offtarget_hits"] == 1
    assert r["blast_weighted_offtarget"] == pytest.approx(1.0)  # exact match → binding=1
    assert r["blast_specificity"] == 0.0


def test_parse_off_target_central_mismatch_downweighted():
    # Central mismatch (pos 15) → binding ≈ 0 → weighted_offtarget ≈ 0.
    # weighted_specificity should be much higher than raw specificity when
    # an off-target near-miss has its mismatch at the probe centre.
    # Need at least one target hit for the ratio to be meaningful.
    target_stitle = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides fragilis"
    off_stitle = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;uncultured bacterium"
    btop = "15GA16"  # mismatch at pos 15 → binding ≈ 0
    blast_out = "\n".join([
        _make_blast_line("kmer_0", target_stitle, 100.0, "32"),
        _make_blast_line("kmer_0", off_stitle, 96.875, btop),
    ])
    results = parse_blast_results(blast_out, [K], "Bacteroides", "genus")
    r = results[K]
    assert r["blast_offtarget_hits"] == 1
    assert r["blast_weighted_offtarget"] == pytest.approx(0.0, abs=0.02)
    # Raw: 1/(1+1) = 0.50; weighted: 1/(1+~0) ≈ 1.0
    assert r["blast_specificity"] == pytest.approx(0.5)
    assert r["blast_weighted_specificity"] > 0.97


def test_parse_off_target_end_mismatch_preserved():
    # mismatch at position 0 (end) → binding score ≈ 0.998
    stitle = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;uncultured bacterium"
    btop = "AG31"  # mismatch at pos 0
    line = _make_blast_line("kmer_0", stitle, 96.875, btop)
    results = parse_blast_results(line, [K], "Bacteroides", "genus")

    r = results[K]
    assert r["blast_weighted_offtarget"] == pytest.approx(0.998, abs=0.01)


def test_parse_uses_qstart_to_locate_mismatch():
    # BTOP "15GA13" starting at query base 4 puts the mismatch at probe pos 18,
    # not 15; bases 0-2 are unaligned. Without the offset this would score as
    # a near-central mismatch.
    stitle = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;uncultured bacterium"
    line = _make_blast_line("kmer_0", stitle, 96.4, "15GA13", length=29) + "\t4\t32"
    results = parse_blast_results(line, [K], "Bacteroides", "genus")

    expected = min(1.0, sum(mismatch_weight(p, 32) for p in [0, 1, 2, 18]))
    assert results[K]["blast_weighted_offtarget"] == pytest.approx(1 - expected, abs=1e-3)


def test_parse_taxonomy_species_with_spaces():
    # Species name "uncultured bacterium" contains a space — must not split on it
    stitle = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;uncultured bacterium"
    line = _make_blast_line("kmer_0", stitle, 100.0, "32")
    results = parse_blast_results(line, [K], "Bacteroides", "genus")
    # Should be on-target, not off-target due to bad space-split
    assert results[K]["blast_target_hits"] == 1


def test_parse_no_blast_output():
    results = parse_blast_results("", [K], "Bacteroides", "genus")
    r = results[K]
    assert r["blast_target_hits"] == 0
    assert r["blast_offtarget_hits"] == 0
    assert r["blast_specificity"] == 1.0  # no hits → treated as fully specific


def test_parse_mixed_hits():
    target_stitle = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides fragilis"
    off_stitle = "Bacteria;Firmicutes;Clostridia;Lachnospirales;Lachnospiraceae;Lachnospiraceae UCG-004;uncultured bacterium"
    blast_out = "\n".join([
        _make_blast_line("kmer_0", target_stitle, 100.0, "32"),
        _make_blast_line("kmer_0", target_stitle, 100.0, "32"),
        _make_blast_line("kmer_0", off_stitle, 100.0, "32"),
    ])
    results = parse_blast_results(blast_out, [K], "Bacteroides", "genus")
    r = results[K]
    assert r["blast_target_hits"] == 2
    assert r["blast_offtarget_hits"] == 1
    assert r["blast_specificity"] == pytest.approx(2 / 3)
    assert r["blast_top_offtarget"] == "Lachnospiraceae UCG-004"


def test_parse_blast_capped_flag():
    stitle = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;B. fragilis"
    lines = [_make_blast_line("kmer_0", stitle, 100.0, "32") for _ in range(1000)]
    blast_out = "\n".join(lines)
    results = parse_blast_results(blast_out, [K], "Bacteroides", "genus", max_target_seqs=1000)
    assert results[K]["blast_capped"] is True


def test_parse_blast_not_capped():
    stitle = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;B. fragilis"
    lines = [_make_blast_line("kmer_0", stitle, 100.0, "32") for _ in range(5)]
    blast_out = "\n".join(lines)
    results = parse_blast_results(blast_out, [K], "Bacteroides", "genus", max_target_seqs=1000)
    assert results[K]["blast_capped"] is False


def test_parse_species_strain_counts_as_target():
    tax = "Bacteria;Campylobacterota;Campylobacteria;Campylobacterales;Helicobacteraceae;Helicobacter;"
    blast_out = "\n".join([
        _make_blast_line("kmer_0", tax + "Helicobacter pylori", 100.0, "32"),
        _make_blast_line("kmer_0", tax + "Helicobacter pylori 26695", 100.0, "32"),
        _make_blast_line("kmer_0", tax + "Helicobacter felis", 100.0, "32"),
    ])
    r = parse_blast_results(blast_out, [K], "Helicobacter pylori", "species")[K]
    assert r["blast_target_hits"] == 2
    assert r["blast_offtarget_hits"] == 1
    assert r["blast_top_offtarget"] == "Helicobacter felis"


def test_parse_species_unnamed_congeners_reported_separately():
    genus_tax = "Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;"
    other_tax = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;"
    blast_out = "\n".join([
        _make_blast_line("kmer_0", genus_tax + "Fusobacterium nucleatum", 100.0, "32"),
        _make_blast_line("kmer_0", genus_tax + "uncultured bacterium", 100.0, "32"),
        _make_blast_line("kmer_0", genus_tax + "Fusobacterium periodonticum", 100.0, "32"),
        _make_blast_line("kmer_0", other_tax + "uncultured bacterium", 100.0, "32"),
    ])
    r = parse_blast_results(blast_out, [K], "Fusobacterium nucleatum", "species",
                            target_genus="Fusobacterium")[K]
    # standard columns unchanged: all 3 non-target hits are off-target
    assert r["blast_offtarget_hits"] == 3
    assert r["blast_specificity"] == pytest.approx(1 / 4)
    # only the unnamed hit in the target genus is set aside
    assert r["blast_unnamed_congener_hits"] == 1
    assert r["blast_specificity_excl_unnamed"] == pytest.approx(1 / 3)
    assert r["blast_weighted_specificity_excl_unnamed"] == pytest.approx(1 / 3)


def test_parse_no_unnamed_keys_without_target_genus():
    stitle = "Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;uncultured bacterium"
    r = parse_blast_results(_make_blast_line("kmer_0", stitle, 100.0, "32"), [K],
                            "Fusobacterium nucleatum", "species")[K]
    assert "blast_unnamed_congener_hits" not in r


# ── blast_batch_size ──────────────────────────────────────────────────────────

def test_blast_batch_size_default_caps_unchanged():
    assert blast_batch_size(1000) == 20
    assert blast_batch_size(50000) == 20


def test_blast_batch_size_shrinks_with_large_caps():
    assert blast_batch_size(200_000) == 10
    assert blast_batch_size(500_000) == 4
    assert blast_batch_size(5_000_000) == 1
