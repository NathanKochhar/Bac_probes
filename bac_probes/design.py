"""Batch probe design: candidate k-mers, exact off-target scoring, optional BLAST,
and selection of a small probe set that reaches a pooled-coverage target.

Used by `bac-probes design`. Everything here is deterministic: whenever two
candidates tie, the tie is broken on the k-mer sequence itself.

Coverage of each candidate is held as a bitmask (a Python int, bit i = target
sequence i), so a probe covering 70,000 sequences costs ~9 KB rather than a
70,000-element set.
"""

from __future__ import annotations

import multiprocessing as mp
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator

import pandas as pd

from .database import TAXONOMY_LEVELS, is_unnamed_species, taxon_matches
from .kmers import (
    extract_kmers,
    gc_content,
    merge_overlapping_kmers,
    rescore_merged_conservation,
    reverse_complement,
)

DEFAULT_SPEC_FLOORS = (1.0, 0.999, 0.995, 0.99, 0.98, 0.95, 0.90)
OVERLAP_WORD = 16  # probes sharing a 16-mer are treated as the same site
SELECTION_LIST_SIZE = 1000  # per ordering, for selection without BLAST


# ── Targets ───────────────────────────────────────────────────────────────────

@dataclass
class Target:
    name: str
    level: str
    exclude_unnamed: bool = False  # species level + --exclude-unnamed-congeners
    seqs: list[str] = field(default_factory=list)
    genera: Counter = field(default_factory=Counter)
    n_ambiguous: int = 0

    @property
    def genus(self) -> str:
        """SILVA genus of the target: the most common genus among its sequences."""
        return self.genera.most_common(1)[0][0] if self.genera else ""

    def classify(self, taxonomy: dict[str, str]) -> str:
        """'target', 'ambiguous' (unnamed same-genus, when excluded) or 'off'."""
        if taxon_matches(taxonomy.get(self.level, ""), self.name, self.level):
            return "target"
        if (
            self.exclude_unnamed
            and taxonomy.get("genus", "") == self.genus
            and is_unnamed_species(taxonomy.get("species", ""))
        ):
            return "ambiguous"
        return "off"


def read_targets(target: str, default_level: str) -> list[tuple[str, str]]:
    """
    Resolve TARGET into (name, level) pairs.

    A path to an existing .csv file is read as a taxon list: it needs a 'taxa'
    column (any case) and may have a 'level' column; rows without a level use
    `default_level`. Anything else is a single taxon name at `default_level`.
    """
    path = Path(target)
    if path.suffix.lower() != ".csv" or not path.is_file():
        return [(target.strip(), default_level)]

    df = pd.read_csv(path, dtype=str)
    cols = {c.strip().lower(): c for c in df.columns}
    if "taxa" not in cols:
        raise ValueError(f"{target}: CSV needs a 'taxa' column (found: {', '.join(df.columns)})")

    out: list[tuple[str, str]] = []
    for _, row in df.iterrows():
        name = row[cols["taxa"]]
        if pd.isna(name) or not str(name).strip():
            continue
        level = row[cols["level"]] if "level" in cols else None
        level = str(level).strip().lower() if level is not None and not pd.isna(level) else default_level
        if level not in TAXONOMY_LEVELS:
            raise ValueError(f"{target}: unknown level '{level}' for '{name}'")
        pair = (str(name).strip(), level)
        if pair not in out:
            out.append(pair)
    if not out:
        raise ValueError(f"{target}: no taxa found")
    return out


# ── Parallel helper ───────────────────────────────────────────────────────────

def run_jobs(
    func: Callable,
    jobs: Iterable,
    threads: int,
    initializer: Callable | None = None,
    initargs: tuple = (),
) -> Iterator:
    """Map `func` over `jobs`, in-process for threads <= 1, else with a process pool."""
    if threads <= 1:
        if initializer:
            initializer(*initargs)
        for job in jobs:
            yield func(job)
        return
    with mp.Pool(threads, initializer=initializer, initargs=initargs) as pool:
        yield from pool.imap_unordered(func, jobs)


# ── Conservation (per target) ─────────────────────────────────────────────────

def conservation_job(args) -> tuple[int, dict[str, int]]:
    """
    Candidate probes for one target: {probe: number of target seqs containing it},
    keeping those present in >= min_conservation of the sequences.
    """
    idx, seqs, k, min_cons, gc_range, max_hp, merge = args
    counts: Counter = Counter()
    for seq in seqs:
        counts.update(extract_kmers(seq, k, gc_range=gc_range, max_homopolymer=max_hp))
    need = max(1.0, min_cons * len(seqs))
    hits = {km: c for km, c in counts.items() if c >= need}
    if merge and hits:
        merged = merge_overlapping_kmers(set(hits), k)
        fracs = rescore_merged_conservation(seqs, merged, min_cons, gc_range=gc_range, max_homopolymer=max_hp)
        hits = {m: round(f * len(seqs)) for m, f in fracs.items()}
    return idx, hits


# ── Exact off-target scan (all targets, one pass) ─────────────────────────────

_SCAN: dict = {}


def init_offtarget_scan(owners: dict[str, list[int]], targets: list[Target], k: int) -> None:
    """Worker initializer: candidate → owning target indices, and the targets (without seqs)."""
    _SCAN["owners"] = owners
    _SCAN["targets"] = targets
    _SCAN["k"] = k
    _SCAN["variable"] = any(len(c) != k for c in owners)


def offtarget_job(chunk: list[tuple[tuple[str, ...], str]]) -> tuple[Counter, Counter]:
    """
    For a chunk of (taxonomy names, sequence) records, count per (target, probe):
      off      – non-target sequences containing the probe
      unnamed  – of those, unnamed sequences in the target's genus (only tallied
                 for targets with exclude_unnamed)
    """
    owners, targets, k = _SCAN["owners"], _SCAN["targets"], _SCAN["k"]
    off: Counter = Counter()
    unnamed: Counter = Counter()
    for names, seq in chunk:
        taxonomy = dict(zip(TAXONOMY_LEVELS, names))
        if _SCAN["variable"]:
            present = [c for c in owners if c in seq]
        else:
            present = owners.keys() & {seq[i:i + k] for i in range(len(seq) - k + 1)}
        for probe in present:
            for tid in owners[probe]:
                cls = targets[tid].classify(taxonomy)
                if cls == "target":
                    continue
                off[(tid, probe)] += 1
                if cls == "ambiguous":
                    unnamed[(tid, probe)] += 1
    return off, unnamed


def chunk_records(records: Iterable[tuple[dict[str, str], str]], size: int = 2000) -> Iterator[list]:
    """Group SILVA records into picklable chunks of (taxonomy names, seq)."""
    chunk: list = []
    for taxonomy, seq in records:
        chunk.append((tuple(taxonomy.get(lv, "") for lv in TAXONOMY_LEVELS), seq))
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


# ── Candidate table ───────────────────────────────────────────────────────────

def candidate_table(
    target: Target,
    hits: dict[str, int],
    off: Counter,
    unnamed: Counter,
    tid: int,
    n_total: int,
) -> pd.DataFrame:
    """One row per candidate with exact coverage and specificity."""
    n_t = len(target.seqs)
    n_bg = n_total - n_t
    n_bg_ex = n_bg - target.n_ambiguous
    rows = []
    for probe_seq, h in hits.items():
        ot = off.get((tid, probe_seq), 0)
        row = {
            "kmer": probe_seq,
            "probe": reverse_complement(probe_seq),
            "length": len(probe_seq),
            "gc_pct": round(gc_content(probe_seq) * 100, 1),
            "conservation_pct": round(100 * h / n_t, 2) if n_t else 0.0,
            "target_seqs": n_t,
            "target_hits": h,
            "exact_offtarget": ot,
            "exact_bg_seqs": n_bg,
            "exact_specificity": round(1 - ot / n_bg, 6) if n_bg else 1.0,
            "exact_precision": round(h / (h + ot), 4) if h + ot else 0.0,
        }
        if target.exclude_unnamed:
            un = unnamed.get((tid, probe_seq), 0)
            ot_ex = ot - un
            row.update({
                "exact_unnamed_congener_hits": un,
                "exact_offtarget_excl_unnamed": ot_ex,
                "exact_specificity_excl_unnamed": round(1 - ot_ex / n_bg_ex, 6) if n_bg_ex else 1.0,
                "exact_precision_excl_unnamed": round(h / (h + ot_ex), 4) if h + ot_ex else 0.0,
            })
        rows.append(row)
    cols = ["kmer", "probe", "length", "gc_pct", "conservation_pct", "target_seqs", "target_hits",
            "exact_offtarget", "exact_bg_seqs", "exact_specificity", "exact_precision"]
    df = pd.DataFrame(rows, columns=cols if not rows else None)
    if df.empty:
        return df
    precision = "exact_precision_excl_unnamed" if target.exclude_unnamed else "exact_precision"
    return df.sort_values([precision, "conservation_pct", "kmer"], ascending=[False, False, True],
                          ignore_index=True)


# ── Coverage bitmasks ─────────────────────────────────────────────────────────

def coverage_masks(seqs: list[str], probes: Iterable[str]) -> dict[str, int]:
    """Bitmask per probe: bit i is set when target sequence i contains the probe."""
    probes = list(dict.fromkeys(probes))
    if not probes:
        return {}
    lengths = {len(p) for p in probes}
    pset = set(probes)
    bufs = {p: bytearray((len(seqs) + 7) // 8) for p in probes}
    for i, seq in enumerate(seqs):
        if len(lengths) == 1:
            k = next(iter(lengths))
            present = pset & {seq[j:j + k] for j in range(len(seq) - k + 1)}
        else:
            present = [p for p in probes if p in seq]
        for p in present:
            bufs[p][i >> 3] |= 1 << (i & 7)
    return {p: int.from_bytes(b, "little") for p, b in bufs.items()}


# ── Selection ─────────────────────────────────────────────────────────────────

def site_words(seq: str) -> set[str]:
    return {seq[i:i + OVERLAP_WORD] for i in range(len(seq) - OVERLAP_WORD + 1)}


def distinct_sites(probes: Iterable[str], limit: int, seen: set[str] | None = None) -> list[str]:
    """First `limit` probes (in the given order) that share no 16-mer with an earlier one."""
    seen = set() if seen is None else seen
    out: list[str] = []
    for p in probes:
        if len(out) >= limit:
            break
        w = site_words(p)
        if w & seen:
            continue
        out.append(p)
        seen |= w
    return out


def _sorted(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    return df.sort_values(cols + ["kmer"], ascending=[False] * len(cols) + [True]).kmer


def blast_pool(df: pd.DataFrame, metric: str, min_floor: float, limit: int) -> list[str]:
    """
    Candidates to BLAST: up to `limit` at distinct sites, half chosen for
    coverage and half for exact specificity, all with `metric` >= min_floor.
    """
    ok = df[df[metric] >= min_floor]
    by_cov = distinct_sites(_sorted(ok, ["conservation_pct"]), limit // 2)
    seen = set().union(*map(site_words, by_cov)) if by_cov else set()
    by_spec = distinct_sites(_sorted(ok, [metric, "conservation_pct"]), limit - len(by_cov), seen)
    return by_cov + by_spec


def selection_pool(df: pd.DataFrame, metric: str, min_floor: float,
                   per_list: int = SELECTION_LIST_SIZE) -> list[str]:
    """Candidates the selection considers without BLAST: the best by specificity and by coverage."""
    ok = df[df[metric] >= min_floor]
    by_spec = list(_sorted(ok, [metric, "conservation_pct"])[:per_list])
    by_cov = list(_sorted(ok, ["conservation_pct"])[:per_list])
    return list(dict.fromkeys(by_spec + by_cov))


def greedy_select(
    df: pd.DataFrame,
    masks: dict[str, int],
    n_target: int,
    metric: str,
    floor: float,
    n_probes: int,
) -> tuple[list[tuple[str, int]], float]:
    """
    Greedy set cover over candidates with `metric` >= floor: each pick adds the
    most not-yet-covered target sequences (ties: higher metric, then k-mer),
    skipping probes at an already-picked site. If coverage stops growing before
    n_probes, the remaining slots are filled with further non-overlapping probes
    (best metric first). Returns [(probe, newly covered seqs)] and pooled fraction.
    """
    pool = df[(df[metric] >= floor) & df.kmer.isin(masks.keys())]
    spec = dict(zip(pool.kmer, pool[metric]))
    remaining = {p: masks[p] for p in pool.kmer}
    covered, seen, picked = 0, set(), []
    while len(picked) < n_probes and remaining:
        best = min(remaining, key=lambda p: (-(remaining[p] & ~covered).bit_count(), -spec[p], p))
        gain = (remaining[best] & ~covered).bit_count()
        if gain == 0:
            break
        mask = remaining.pop(best)
        if site_words(best) & seen:
            continue
        picked.append((best, gain))
        covered |= mask
        seen |= site_words(best)
    for p in _sorted(pool, [metric, "conservation_pct"]):
        if len(picked) >= n_probes:
            break
        if p in remaining and not site_words(p) & seen:
            picked.append((p, 0))
            seen |= site_words(p)
    return picked, (covered.bit_count() / n_target if n_target else 0.0)


def select_probes(
    df: pd.DataFrame,
    masks: dict[str, int],
    n_target: int,
    metric: str,
    floors: Iterable[float],
    n_probes: int,
    pool_target: float,
) -> tuple[float | None, list[tuple[str, int]], float]:
    """
    Try specificity floors from strict to loose; return the first whose probe
    set reaches pool_target, else the floor giving the highest pooled coverage.
    """
    best: tuple[float | None, list, float] = (None, [], 0.0)
    for floor in sorted(set(floors), reverse=True):
        picked, pooled = greedy_select(df, masks, n_target, metric, floor, n_probes)
        if best[0] is None or pooled > best[2] + 1e-12:
            best = (floor, picked, pooled)
        if pooled >= pool_target:
            return floor, picked, pooled
    return best


def selected_table(df: pd.DataFrame, picked: list[tuple[str, int]], n_target: int) -> pd.DataFrame:
    """Rows of `df` for the picked probes, in pick order, with coverage bookkeeping."""
    if not picked:
        return df.iloc[0:0].assign(probe_rank=pd.Series(dtype=int))
    sel = df.set_index("kmer").loc[[p for p, _ in picked]].reset_index()
    sel.insert(0, "probe_rank", range(1, len(sel) + 1))
    sel["new_seqs_covered"] = [n for _, n in picked]
    sel["cumulative_coverage_pct"] = (100 * sel.new_seqs_covered.cumsum() / n_target).round(2)
    return sel


def safe_filename(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name.strip())
