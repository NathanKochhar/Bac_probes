"""Tests for bac_probes.design and the `design` CLI command (exact scoring, no BLAST)."""

import random
from collections import Counter

import pandas as pd
import pytest
from click.testing import CliRunner

from bac_probes.cli import design
from bac_probes.design import (
    Target,
    blast_pool,
    conservation_job,
    coverage_masks,
    distinct_sites,
    SpacingChecker,
    greedy_select,
    init_offtarget_scan,
    offtarget_job,
    read_targets,
    select_probes,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _rand_seq(rng: random.Random, n: int) -> str:
    """Random sequence with GC ~50% and no homopolymer runs of 4+."""
    out = []
    while len(out) < n:
        b = rng.choice("ACGT")
        if len(out) >= 3 and out[-1] == out[-2] == out[-3] == b:
            continue
        out.append(b)
    return "".join(out)


RNG = random.Random(7)
CORE1 = _rand_seq(RNG, 60)
CORE2 = _rand_seq(RNG, 60)
SHARED = _rand_seq(RNG, 60)  # also present outside the target genus

ALPHA = "Bacteria;P1;C1;O1;F1;Alphagenus;"
BETA = "Bacteria;P2;C2;O2;F2;Betagenus;"


def _write_fasta(path, records):
    path.write_text("".join(f">{acc} {tax}\n{seq}\n" for acc, tax, seq in records))
    return str(path)


@pytest.fixture()
def silva(tmp_path):
    """
    Alphagenus: a1-a2 carry CORE1, a3-a4 carry CORE2, all four carry SHARED.
    Betagenus:  b1 carries SHARED (so SHARED k-mers have an off-target).
    """
    rng = random.Random(11)
    recs = [
        ("a1", ALPHA + "Alphagenus alpha", _rand_seq(rng, 40) + CORE1 + SHARED + _rand_seq(rng, 40)),
        ("a2", ALPHA + "Alphagenus alpha 2", _rand_seq(rng, 40) + CORE1 + SHARED + _rand_seq(rng, 40)),
        ("a3", ALPHA + "Alphagenus gamma", _rand_seq(rng, 40) + CORE2 + SHARED + _rand_seq(rng, 40)),
        ("a4", ALPHA + "Alphagenus gamma", _rand_seq(rng, 40) + CORE2 + SHARED + _rand_seq(rng, 40)),
        ("b1", BETA + "Betagenus beta", _rand_seq(rng, 40) + SHARED + _rand_seq(rng, 40)),
        ("b2", BETA + "Betagenus beta", _rand_seq(rng, 160)),
    ]
    return _write_fasta(tmp_path / "silva.fasta", recs)


def _run(args):
    res = CliRunner().invoke(design, args)
    assert res.exit_code == 0, res.output
    return res


# ── read_targets ──────────────────────────────────────────────────────────────

def test_read_targets_single_name():
    assert read_targets("Fusobacterium", "genus") == [("Fusobacterium", "genus")]


def test_read_targets_csv(tmp_path):
    p = tmp_path / "taxa.csv"
    p.write_text("Taxa,Level\nStaphylococcus,Genus\nEscherichia coli,Species\nStaphylococcus,Genus\n")
    assert read_targets(str(p), "genus") == [("Staphylococcus", "genus"), ("Escherichia coli", "species")]


def test_read_targets_csv_without_level_uses_default(tmp_path):
    p = tmp_path / "taxa.csv"
    p.write_text("taxa\nFirmicutes\n")
    assert read_targets(str(p), "phylum") == [("Firmicutes", "phylum")]


def test_read_targets_csv_needs_taxa_column(tmp_path):
    p = tmp_path / "taxa.csv"
    p.write_text("name,level\nFirmicutes,phylum\n")
    with pytest.raises(ValueError, match="taxa"):
        read_targets(str(p), "genus")


def test_read_targets_csv_bad_level(tmp_path):
    p = tmp_path / "taxa.csv"
    p.write_text("taxa,level\nFirmicutes,kingdom\n")
    with pytest.raises(ValueError, match="kingdom"):
        read_targets(str(p), "genus")


# ── Target.classify ───────────────────────────────────────────────────────────

def _tax(genus, species):
    return {"domain": "Bacteria", "genus": genus, "species": species}


def test_classify_species_with_unnamed_exclusion():
    t = Target("Alphagenus alpha", "species", exclude_unnamed=True, genera=Counter({"Alphagenus": 3}))
    assert t.classify(_tax("Alphagenus", "Alphagenus alpha 2")) == "target"
    assert t.classify(_tax("Alphagenus", "uncultured bacterium")) == "ambiguous"
    assert t.classify(_tax("Alphagenus", "Alphagenus gamma")) == "off"
    assert t.classify(_tax("Betagenus", "uncultured bacterium")) == "off"


def test_classify_without_exclusion():
    t = Target("Alphagenus alpha", "species", genera=Counter({"Alphagenus": 3}))
    assert t.classify(_tax("Alphagenus", "uncultured bacterium")) == "off"


# ── conservation / off-target / coverage ──────────────────────────────────────

def test_conservation_job_threshold():
    seqs = [CORE1 + CORE2, CORE1, CORE1]
    _, hits = conservation_job((0, seqs, 32, 0.5, (0.0, 1.0), 99, False))
    k = CORE1[:32]
    assert hits[k] == 3
    assert CORE2[:32] not in hits  # in 1 of 3 seqs < 50%


def test_offtarget_job_counts_and_unnamed():
    probe = CORE1[:32]
    t = Target("Alphagenus alpha", "species", exclude_unnamed=True, genera=Counter({"Alphagenus": 2}))
    init_offtarget_scan({probe: [0]}, [t], 32)
    chunk = [
        (("Bacteria", "", "", "", "", "Alphagenus", "Alphagenus alpha"), CORE1),        # target
        (("Bacteria", "", "", "", "", "Alphagenus", "uncultured bacterium"), CORE1),    # ambiguous
        (("Bacteria", "", "", "", "", "Betagenus", "Betagenus beta"), "AA" + CORE1),    # off
        (("Bacteria", "", "", "", "", "Betagenus", "Betagenus beta"), CORE2),           # no hit
    ]
    off, unnamed = offtarget_job(chunk)
    assert off[(0, probe)] == 2
    assert unnamed[(0, probe)] == 1


def test_coverage_masks_fixed_and_variable_length():
    seqs = [CORE1, CORE2, CORE1 + CORE2]
    m = coverage_masks(seqs, [CORE1[:32], CORE2[:32]])
    assert m[CORE1[:32]] == 0b101
    assert m[CORE2[:32]] == 0b110
    m2 = coverage_masks(seqs, [CORE1[:40], CORE2[:32]])  # mixed lengths → substring search
    assert m2[CORE1[:40]] == 0b101


# ── Selection ─────────────────────────────────────────────────────────────────

def test_distinct_sites_skips_shared_16mer():
    a, b, c = CORE1[0:32], CORE1[10:42], CORE2[0:32]
    assert distinct_sites([a, b, c], 5) == [a, c]


def _cands(rows):
    return pd.DataFrame(rows, columns=["kmer", "conservation_pct", "exact_precision"])


def test_greedy_picks_complementary_probes():
    p1, p2, p3 = CORE1[:32], CORE2[:32], SHARED[:32]
    df = _cands([(p1, 50, 1.0), (p2, 50, 1.0), (p3, 75, 1.0)])
    masks = {p1: 0b0011, p2: 0b1100, p3: 0b0111}
    picked, pooled = greedy_select(df, masks, 4, "exact_precision", 1.0, 10)
    assert [p for p, _ in picked][:2] == [p3, p2]  # p3 covers most, then p2 adds the rest
    assert pooled == 1.0


def test_greedy_top_up_fills_remaining_slots():
    p1, p2 = CORE1[:32], CORE2[:32]
    df = _cands([(p1, 100, 1.0), (p2, 100, 0.99)])
    masks = {p1: 0b11, p2: 0b11}
    picked, pooled = greedy_select(df, masks, 2, "exact_precision", 0.9, 5)
    assert picked == [(p1, 2), (p2, 0)]  # p2 adds nothing new but tops up
    assert pooled == 1.0


def test_greedy_tie_breaks_on_kmer():
    pa, pb = sorted([CORE1[:32], CORE2[:32]])
    df = _cands([(pb, 50, 1.0), (pa, 50, 1.0)])
    masks = {pa: 0b01, pb: 0b10}
    picked, _ = greedy_select(df, masks, 2, "exact_precision", 1.0, 1)
    assert picked == [(pa, 1)]


def test_select_probes_uses_strictest_floor_reaching_target():
    strict, loose = CORE1[:32], CORE2[:32]
    df = _cands([(strict, 25, 1.0), (loose, 100, 0.95)])
    masks = {strict: 0b0001, loose: 0b1111}
    floor, picked, pooled = select_probes(df, masks, 4, "exact_precision", [1.0, 0.95], 10, 0.5)
    assert floor == 0.95 and pooled == 1.0
    floor, picked, pooled = select_probes(df, masks, 4, "exact_precision", [1.0, 0.95], 10, 0.2)
    assert floor == 1.0 and [p for p, _ in picked] == [strict]


def test_select_probes_falls_back_to_best_coverage():
    p = CORE1[:32]
    df = _cands([(p, 25, 0.5)])
    floor, picked, pooled = select_probes(df, {p: 0b1}, 4, "exact_precision", [1.0, 0.9], 10, 0.5)
    assert picked == [] and pooled == 0.0


def test_blast_pool_mixes_coverage_and_specificity():
    rows = [(CORE1[i:i + 32], 90 - i, 0.95) for i in range(0, 20, 20)]
    rows += [(CORE2[:32], 10, 1.0), (SHARED[:32], 5, 0.5)]
    df = _cands(rows)
    pool = blast_pool(df, "exact_precision", 0.9, 2)
    assert pool == [CORE1[:32], CORE2[:32]]  # one by coverage, one by specificity; SHARED below floor


# ── CLI ───────────────────────────────────────────────────────────────────────

def test_design_single_taxon(silva, tmp_path):
    out = tmp_path / "out"
    _run([silva, "Alphagenus", "--level", "genus", "--threads", "1", "-o", str(out)])
    summary = pd.read_csv(out / "design_summary.csv")
    row = summary.iloc[0]
    assert row.target_seqs == 4
    assert row.spec_floor_used == 1.0  # CORE1 + CORE2 probes are exact-specific
    assert row.pooled_coverage_pct == 100.0
    assert bool(row.goal_reached)
    sel = pd.read_csv(out / "selected" / "Alphagenus_selected_probes.csv")
    assert (sel.exact_offtarget == 0).all()
    assert sel.kmer.str.len().eq(32).all()
    assert (out / "candidates" / "Alphagenus_candidates.tsv").exists()
    cmd = (out / "command.txt").read_text()
    assert "pool_target = 0.5" in cmd and "n_probes = 10" in cmd


def test_design_candidates_include_offtargets(silva, tmp_path):
    out = tmp_path / "out"
    _run([silva, "Alphagenus", "--threads", "1", "--no-select", "-o", str(out)])
    cands = pd.read_csv(out / "candidates" / "Alphagenus_candidates.tsv", sep="\t")
    shared = cands[cands.kmer == SHARED[:32]].iloc[0]
    assert shared.target_hits == 4 and shared.exact_offtarget == 1
    assert shared.exact_precision == pytest.approx(0.8)
    assert not (out / "design_summary.csv").exists()


def test_design_pool_target_and_n_probes(silva, tmp_path):
    out = tmp_path / "out"
    _run([silva, "Alphagenus", "--threads", "1", "--n-probes", "1", "--spec-floors", "1.0", "-o", str(out)])
    row = pd.read_csv(out / "design_summary.csv").iloc[0]
    assert row.probes_selected == 1
    assert row.pooled_coverage_pct == 50.0  # one specific probe covers one of the two cores


def test_design_csv_of_taxa_multiprocess(silva, tmp_path):
    taxa = tmp_path / "taxa.csv"
    taxa.write_text("Taxa,Level\nAlphagenus,Genus\nBetagenus,Genus\nNotagenus,Genus\n")
    out = tmp_path / "out"
    res = _run([silva, str(taxa), "--threads", "2", "-o", str(out)])
    assert "Notagenus" in res.output  # warned and skipped
    summary = pd.read_csv(out / "design_summary.csv")
    assert list(summary.taxa) == ["Alphagenus", "Betagenus"]
    allsel = pd.read_csv(out / "all_selected_probes.csv")
    assert set(allsel.taxa) == {"Alphagenus", "Betagenus"}


def test_design_species_exclude_unnamed(tmp_path):
    rng = random.Random(3)
    recs = [
        ("a1", ALPHA + "Alphagenus alpha", _rand_seq(rng, 30) + CORE1 + _rand_seq(rng, 30)),
        ("a2", ALPHA + "Alphagenus alpha str. X", _rand_seq(rng, 30) + CORE1 + _rand_seq(rng, 30)),
        ("u1", ALPHA + "uncultured bacterium", _rand_seq(rng, 30) + CORE1 + _rand_seq(rng, 30)),
        ("g1", ALPHA + "Alphagenus gamma", _rand_seq(rng, 120)),
    ]
    silva = _write_fasta(tmp_path / "s.fasta", recs)
    out = tmp_path / "out"
    _run([silva, "Alphagenus alpha", "--level", "species", "--exclude-unnamed-congeners",
          "--threads", "1", "-o", str(out)])
    row = pd.read_csv(out / "design_summary.csv").iloc[0]
    assert row.target_seqs == 2                      # strain counts as the species
    assert row.silva_genus == "Alphagenus"
    assert row.ambiguous_seqs_excluded == 1
    assert row.spec_metric == "exact_precision_excl_unnamed"
    sel = pd.read_csv(out / "selected" / "Alphagenus_alpha_selected_probes.csv")
    core = sel[sel.kmer.isin([CORE1[i:i + 32] for i in range(29)])]
    assert not core.empty
    assert (core.exact_offtarget == 1).all() and (core.exact_precision_excl_unnamed == 1.0).all()


def test_design_target_fasta(silva, tmp_path):
    fa = tmp_path / "targets.fasta"
    fa.write_text(f">t1\n{CORE1}\n>t2\n{CORE1}\n")
    out = tmp_path / "out"
    _run([silva, "Alphagenus", "--target-fasta", str(fa), "--threads", "1", "-o", str(out)])
    row = pd.read_csv(out / "design_summary.csv").iloc[0]
    assert row.target_seqs == 2 and row.pooled_coverage_pct == 100.0


def test_design_unknown_taxon_fails(silva, tmp_path):
    res = CliRunner().invoke(design, [silva, "Nosuchgenus", "--threads", "1", "-o", str(tmp_path / "o")])
    assert res.exit_code != 0


def test_design_merge_overlapping(silva, tmp_path):
    out = tmp_path / "out"
    _run([silva, "Alphagenus", "--merge-overlapping", "--threads", "1", "-o", str(out)])
    sel = pd.read_csv(out / "selected" / "Alphagenus_selected_probes.csv")
    assert (sel.length > 32).any()                  # merged probes are longer than k
    # the probes that add coverage sit in the two genus-specific cores
    gainers = sel[sel.new_seqs_covered > 0].kmer
    assert len(gainers) == 2
    assert all(any(s in seg or seg in s for seg in (CORE1, CORE2)) for s in gainers)
    assert sel.cumulative_coverage_pct.iloc[-1] == 100.0


# ── --min-gap spacing ─────────────────────────────────────────────────────────

A, B, C = CORE1[:32], CORE2[:32], SHARED[:32]
F5, F20, F30 = "ACGTA", "ACGTA" * 4, "ACGTA" * 6


def test_spacing_min_spacing_gap_overlap_and_absent():
    seqs = [A + F5 + B, A[:10] + B, C]
    sc = SpacingChecker(seqs, 10)
    assert sc.min_spacing(A, B) == 5                       # seq 0: 5 bp apart
    assert sc.min_spacing(A, C) is None                    # never on the same sequence
    overlap = SpacingChecker([CORE1], 10)
    assert overlap.min_spacing(CORE1[0:32], CORE1[20:52]) == -12   # 12 bp overlap


def test_spacing_conflicts_threshold():
    sc = SpacingChecker([A + F5 + B], 10)
    assert sc.conflicts(B, [A])
    assert not SpacingChecker([A + F5 + B], 5).conflicts(B, [A])   # gap 5 >= 5
    assert not SpacingChecker([A + F5 + B], 0).conflicts(B, [A])


def test_spacing_allows_variants_on_different_sequences():
    sc = SpacingChecker([A, CORE1[2:34]], 10)                # overlapping site, different seqs
    assert not sc.conflicts(CORE1[2:34], [A])


def _spacing_case(min_gap):
    seqs = [A + F5 + B + F20 + C, A + F30 + C, B, A]
    df = _cands([(A, 75, 1.0), (B, 50, 1.0), (C, 50, 1.0)])
    masks = coverage_masks(seqs, [A, B, C])
    return greedy_select(df, masks, len(seqs), "exact_precision", 1.0, 3,
                         SpacingChecker(seqs, min_gap))


def test_greedy_respects_min_gap():
    picked, _ = _spacing_case(10)
    probes = [p for p, _ in picked]
    assert B not in probes                 # only 5 bp from A on seq 0
    assert probes == [A, C]                # C is >= 30 bp from A everywhere


def test_greedy_min_gap_zero_allows_close_probes():
    picked, _ = _spacing_case(0)
    assert [p for p, _ in picked][:2] == [A, B]


def test_design_selected_probes_respect_min_gap(silva, tmp_path):
    out = tmp_path / "out"
    _run([silva, "Alphagenus", "--threads", "1", "-o", str(out)])
    sel = pd.read_csv(out / "selected" / "Alphagenus_selected_probes.csv")
    assert "min_spacing_bp" not in sel.columns
    seqs = [rec.splitlines()[1] for rec in open(silva).read().split(">")[1:]
            if "Alphagenus" in rec.splitlines()[0]]
    sc = SpacingChecker(seqs, 10)
    probes = list(sel.kmer)
    for i, p in enumerate(probes):
        assert not sc.conflicts(p, probes[:i] + probes[i + 1:])
    assert "min_gap = 10" in (out / "command.txt").read_text()
