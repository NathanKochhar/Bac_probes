"""Tests for the validate CLI command (exact matching, no BLAST)."""

import pandas as pd
import pytest
from click.testing import CliRunner

from bac_probes.cli import validate

GENUS = "Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;"
OTHER = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;"
PROBE = "ACGTACGTACGTACGTACGTACGTACGTACGT"
HIT = "A" * 10 + PROBE + "T" * 10
MISS = "C" * 52


@pytest.fixture()
def silva(tmp_path):
    records = [
        ("a1", GENUS + "Fusobacterium nucleatum", HIT),
        ("a2", GENUS + "Fusobacterium nucleatum subsp. animalis", HIT),   # strain → target
        ("a3", GENUS + "Fusobacterium nucleatum", MISS),
        ("a4", GENUS + "uncultured bacterium", HIT),                      # unnamed congener
        ("a5", GENUS + "Fusobacterium periodonticum", HIT),               # named other species
        ("a6", OTHER + "uncultured bacterium", HIT),                      # unnamed, other genus
        ("a7", OTHER + "Lactobacillus acidophilus", MISS),
    ]
    p = tmp_path / "silva.fasta"
    p.write_text("".join(f">{a} {t}\n{s}\n" for a, t, s in records))
    return str(p)


@pytest.fixture()
def probes_csv(tmp_path):
    p = tmp_path / "probes.csv"
    p.write_text(f"Fusobacterium nucleatum,species,{PROBE}\n")
    return str(p)


def _run(silva, probes_csv, tmp_path, *extra):
    out = tmp_path / "out.tsv"
    res = CliRunner().invoke(
        validate, [silva, probes_csv, "--no-header", "--output", str(out), *extra]
    )
    assert res.exit_code == 0, res.output
    return pd.read_csv(out, sep="\t").iloc[0]


def test_validate_counts_strains_as_target(silva, probes_csv, tmp_path):
    r = _run(silva, probes_csv, tmp_path)
    assert r.target_seqs == 3
    assert r.coverage_count == 2
    assert r.exact_offtarget == 3          # unnamed congener, F. periodonticum, other-genus unnamed
    assert "target_genus" not in r.index   # extra columns only with the flag


def test_validate_exclude_unnamed_congeners(silva, probes_csv, tmp_path):
    r = _run(silva, probes_csv, tmp_path, "--exclude-unnamed-congeners")
    assert r.target_genus == "Fusobacterium"
    assert r.exact_offtarget == 3
    assert r.exact_unnamed_congener_hits == 1
    assert r.exact_offtarget_excl_unnamed == 2
    # background without the unnamed congener: 4 → 3 seqs, 2 of them hit
    assert r.exact_specificity_excl_unnamed == pytest.approx(1 - 2 / 3, abs=1e-6)


# ── --breakdown ───────────────────────────────────────────────────────────────

@pytest.fixture()
def genus_probes_csv(tmp_path):
    p = tmp_path / "genus_probes.csv"
    p.write_text(f"Fusobacterium,genus,{PROBE}\n")
    return str(p)


def test_validate_breakdown_table(silva, genus_probes_csv, tmp_path):
    out = tmp_path / "res.tsv"
    res = CliRunner().invoke(validate, [
        silva, genus_probes_csv, "--no-header", "--breakdown", "species",
        "--min-seqs", "1", "--output", str(out),
    ])
    assert res.exit_code == 0, res.output
    bd = pd.read_csv(tmp_path / "res_breakdown.tsv", sep="\t")
    assert list(bd.columns) == ["taxa", "level", "subtaxon", "total_seqs", "k1", "pool_coverage"]
    rows = bd.set_index("subtaxon")
    assert rows.loc["ALL", "total_seqs"] == 5
    assert rows.loc["ALL", "pool_coverage"] == "4/5 (80.0%)"
    assert rows.loc["Fusobacterium nucleatum", "k1"] == "1/2 (50%)"
    assert rows.loc["Fusobacterium periodonticum", "k1"] == "1/1 (100%)"
    assert "k1 = " + PROBE in res.output


def test_validate_breakdown_min_seqs(silva, genus_probes_csv, tmp_path):
    res = CliRunner().invoke(validate, [
        silva, genus_probes_csv, "--no-header", "--breakdown", "species", "--min-seqs", "2",
    ])
    assert res.exit_code == 0, res.output
    assert "Fusobacterium nucleatum " in res.output   # 2 seqs → shown
    assert "periodonticum" not in res.output.split("Coverage by species")[1]  # 1 seq → hidden


def test_validate_breakdown_multiple_probes_and_taxa(silva, tmp_path):
    other = "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"
    p = tmp_path / "multi.csv"
    p.write_text(f"Fusobacterium,genus,{PROBE}\nFusobacterium,genus,{other}\nLactobacillus,genus,{PROBE}\n")
    out = tmp_path / "multi.tsv"
    res = CliRunner().invoke(validate, [silva, str(p), "--no-header", "--breakdown", "species",
                                        "--min-seqs", "1", "-o", str(out)])
    assert res.exit_code == 0, res.output
    bd = pd.read_csv(tmp_path / "multi_breakdown.tsv", sep="\t", keep_default_na=False)
    fuso_all = bd[(bd.taxa == "Fusobacterium") & (bd.subtaxon == "ALL")].iloc[0]
    assert fuso_all.k2 == "0/5 (0.0%)"
    lacto = bd[bd.taxa == "Lactobacillus"]
    assert set(lacto.subtaxon) == {"ALL", "uncultured bacterium", "Lactobacillus acidophilus"}
    assert (lacto.k2 == "").all()  # Lactobacillus has a single probe


def test_validate_without_breakdown_writes_no_breakdown_file(silva, genus_probes_csv, tmp_path):
    out = tmp_path / "res.tsv"
    res = CliRunner().invoke(validate, [silva, genus_probes_csv, "--no-header", "--output", str(out)])
    assert res.exit_code == 0, res.output
    assert not (tmp_path / "res_breakdown.tsv").exists()
