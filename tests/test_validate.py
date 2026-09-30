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
