"""Unit tests for bac_probes.database."""

import pytest
from bac_probes.database import parse_taxonomy, TAXONOMY_LEVELS


# ── parse_taxonomy ────────────────────────────────────────────────────────────

def test_parse_full_taxonomy():
    s = "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides fragilis"
    t = parse_taxonomy(s)
    assert t["domain"] == "Bacteria"
    assert t["phylum"] == "Bacteroidota"
    assert t["class"] == "Bacteroidia"
    assert t["order"] == "Bacteroidales"
    assert t["family"] == "Bacteroidaceae"
    assert t["genus"] == "Bacteroides"
    assert t["species"] == "Bacteroides fragilis"

def test_parse_taxonomy_with_spaces_in_species():
    # Species names often contain spaces ("uncultured bacterium")
    s = "Bacteria;Firmicutes;Clostridia;Lachnospirales;Lachnospiraceae;Lachnospiraceae UCG-004;uncultured bacterium"
    t = parse_taxonomy(s)
    assert t["genus"] == "Lachnospiraceae UCG-004"
    assert t["species"] == "uncultured bacterium"

def test_parse_taxonomy_incomplete():
    s = "Bacteria;Fusobacteriota"
    t = parse_taxonomy(s)
    assert t["domain"] == "Bacteria"
    assert t["phylum"] == "Fusobacteriota"
    assert t["class"] == ""
    assert t["genus"] == ""

def test_parse_taxonomy_empty_string():
    t = parse_taxonomy("")
    for level in TAXONOMY_LEVELS:
        assert t[level] == ""

def test_parse_taxonomy_all_levels_present():
    t = parse_taxonomy("A;B;C;D;E;F;G")
    assert t["domain"] == "A"
    assert t["species"] == "G"

def test_parse_taxonomy_strips_whitespace():
    s = "Bacteria; Firmicutes ; Bacilli"
    t = parse_taxonomy(s)
    assert t["phylum"] == "Firmicutes"
    assert t["class"] == "Bacilli"

def test_parse_taxonomy_returns_all_levels():
    t = parse_taxonomy("Bacteria;Firmicutes")
    assert set(t.keys()) == set(TAXONOMY_LEVELS)

def test_taxonomy_levels_order():
    assert TAXONOMY_LEVELS[0] == "domain"
    assert TAXONOMY_LEVELS[-1] == "species"
    assert "genus" in TAXONOMY_LEVELS
    assert "phylum" in TAXONOMY_LEVELS
