"""Unit tests for bac_probes.database."""

import pytest
from bac_probes.database import is_unnamed_species, parse_taxonomy, taxon_matches, TAXONOMY_LEVELS


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


# ── taxon_matches ─────────────────────────────────────────────────────────────

def test_taxon_matches_exact_case_insensitive():
    assert taxon_matches("Helicobacter pylori", "helicobacter pylori", "species")


def test_taxon_matches_species_strain():
    assert taxon_matches("Helicobacter pylori 26695", "Helicobacter pylori", "species")


def test_taxon_matches_species_subspecies():
    assert taxon_matches("Campylobacter jejuni subsp. doylei", "Campylobacter jejuni", "species")


def test_taxon_matches_species_rejects_other_species():
    assert not taxon_matches("Helicobacter felis", "Helicobacter pylori", "species")


def test_taxon_matches_species_needs_word_boundary():
    # "Escherichia colix" is not a strain of "Escherichia coli"
    assert not taxon_matches("Escherichia colix", "Escherichia coli", "species")


def test_taxon_matches_subspecies_target_is_exactish():
    target = "Enterobacter hormaechei subsp. oharae"
    assert taxon_matches("Enterobacter hormaechei subsp. oharae ECR091", target, "species")
    assert not taxon_matches("Enterobacter hormaechei subsp. hoffmannii", target, "species")


def test_taxon_matches_genus_is_exact_only():
    assert taxon_matches("Clostridium", "Clostridium", "genus")
    assert not taxon_matches("Clostridium sensu stricto 1", "Clostridium", "genus")


def test_taxon_matches_empty():
    assert not taxon_matches("", "Helicobacter pylori", "species")
    assert not taxon_matches("Helicobacter pylori", "", "species")


# ── is_unnamed_species ────────────────────────────────────────────────────────

@pytest.mark.parametrize("name", [
    "uncultured bacterium",
    "uncultured Helicobacter sp.",
    "Fusobacterium sp.",
    "Fusobacterium sp. oral taxon 203",
    "unidentified",
    "metagenome",
    "",
])
def test_is_unnamed_species_true(name):
    assert is_unnamed_species(name)


@pytest.mark.parametrize("name", [
    "Helicobacter pylori",
    "Helicobacter pylori 26695",
    "Campylobacter jejuni subsp. doylei",
    "Fusobacterium periodonticum",
])
def test_is_unnamed_species_false(name):
    assert not is_unnamed_species(name)
