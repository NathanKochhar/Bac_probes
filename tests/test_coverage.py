"""Unit tests for the coverage CLI command."""

import gzip
import os
import tempfile

import pytest
from click.testing import CliRunner

from bac_probes.cli import coverage


def _write_fasta(path: str, records: list[tuple[str, str, str]]) -> None:
    """Write a tiny SILVA-style FASTA (header includes taxonomy after first space)."""
    with open(path, "w") as fh:
        for acc, taxonomy, seq in records:
            fh.write(f">{acc} {taxonomy}\n{seq}\n")


FUSO_TAX = "Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium nucleatum"
OTHER_TAX = "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus acidophilus"
KMER = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32-mer

FUSO_SEQ_HIT = "A" * 10 + KMER + "T" * 10
FUSO_SEQ_MISS = "C" * 52
OTHER_SEQ = "G" * 52


@pytest.fixture()
def tiny_fasta(tmp_path):
    p = tmp_path / "tiny.fasta"
    _write_fasta(str(p), [
        ("acc1.1.52", FUSO_TAX, FUSO_SEQ_HIT),
        ("acc2.1.52", FUSO_TAX, FUSO_SEQ_MISS),
        ("acc3.1.52", FUSO_TAX, FUSO_SEQ_HIT),
        ("acc4.1.52", OTHER_TAX, OTHER_SEQ),
    ])
    return str(p)


BASE_ARGS = ["--level", "genus", "--breakdown", "species", "--kmer", KMER, "--min-seqs", "1"]


def test_coverage_basic(tiny_fasta):
    runner = CliRunner()
    result = runner.invoke(coverage, [tiny_fasta, "Fusobacterium"] + BASE_ARGS)
    assert result.exit_code == 0, result.output
    # Should show 2/3 hits for Fusobacterium nucleatum
    assert "2/3" in result.output
    assert "Fusobacterium nucleatum" in result.output


def test_coverage_pool_column(tiny_fasta):
    runner = CliRunner()
    result = runner.invoke(coverage, [tiny_fasta, "Fusobacterium"] + BASE_ARGS)
    assert result.exit_code == 0
    assert "pool_coverage" in result.output


def test_coverage_all_summary_row(tiny_fasta):
    runner = CliRunner()
    result = runner.invoke(coverage, [tiny_fasta, "Fusobacterium"] + BASE_ARGS)
    assert result.exit_code == 0
    assert "ALL" in result.output


def test_coverage_output_tsv(tiny_fasta, tmp_path):
    out = str(tmp_path / "cov.tsv")
    runner = CliRunner()
    result = runner.invoke(coverage, [
        tiny_fasta, "Fusobacterium",
        "--level", "genus",
        "--breakdown", "species",
        "--kmer", KMER,
        "--min-seqs", "1",
        "--output", out,
    ])
    assert result.exit_code == 0
    assert os.path.exists(out)
    import pandas as pd
    df = pd.read_csv(out, sep="\t")
    assert "subtaxon" in df.columns
    assert "pool_coverage" in df.columns


def test_coverage_min_seqs_filter(tiny_fasta):
    runner = CliRunner()
    result = runner.invoke(coverage, [
        tiny_fasta, "Fusobacterium",
        "--level", "genus",
        "--breakdown", "species",
        "--kmer", KMER,
        "--min-seqs", "100",  # filter out everything
    ])
    assert result.exit_code == 0
    # Only the ALL row should remain (not filtered by min-seqs)
    assert "ALL" in result.output


def test_coverage_no_kmers_error(tiny_fasta):
    runner = CliRunner()
    result = runner.invoke(coverage, [tiny_fasta, "Fusobacterium"])
    assert result.exit_code != 0


def test_coverage_unknown_taxon(tiny_fasta):
    runner = CliRunner()
    result = runner.invoke(coverage, [
        tiny_fasta, "Streptococcus",
        "--level", "genus",
        "--kmer", KMER,
    ])
    assert result.exit_code != 0


def test_coverage_from_tsv(tiny_fasta, tmp_path):
    import pandas as pd
    probes_tsv = str(tmp_path / "probes.tsv")
    pd.DataFrame([
        {"kmer": KMER, "exact_specificity": 1.0, "blast_weighted_specificity": 0.99},
        {"kmer": "C" * 32, "exact_specificity": 0.9, "blast_weighted_specificity": 0.8},
    ]).to_csv(probes_tsv, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(coverage, [
        tiny_fasta, "Fusobacterium",
        "--level", "genus",
        "--from-tsv", probes_tsv,
        "--top-n", "1",
    ])
    assert result.exit_code == 0
    assert KMER in result.output
