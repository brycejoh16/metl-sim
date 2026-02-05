#!/usr/bin/env python3
"""
mutation_subset_picker.py

Updates:
- summary.json now includes per-k progress:
    completed / total / remaining / percent
  both for all ks and for ks within --k_range.

Run tests:
  python -m unittest -v mutation_subset_picker.py

Arg-file support:
  python mutation_subset_picker.py @args.txt
"""

from __future__ import annotations

import argparse
import json
import random
import unittest
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None  # type: ignore


# =========================
# Core logic
# =========================

@dataclass(frozen=True, order=True)
class MutCombo:
    pdb: str
    muts: Tuple[str, ...]  # sorted tuple

    def to_line(self) -> str:
        return f"{self.pdb}\t{','.join(self.muts)}"


def parse_k_range_arg(s: Optional[str]) -> Tuple[str, int, int]:
    """
    Parse --k_range.

    Returns:
      ("all", 1, 1) if s is None or empty/whitespace
      ("range", a, b) if s is "a-b"
      ("single_k", k, k) if s is integer like "2"
    """
    if s is None:
        return ("all", 1, 1)
    s = s.strip()
    if s == "":
        return ("all", 1, 1)

    if "-" in s:
        parts = s.split("-")
        if len(parts) != 2:
            raise ValueError(f"Invalid --k_range: {s!r}. Use like 2-3.")
        a = int(parts[0])
        b = int(parts[1])
        if a <= 0 or b <= 0 or b < a:
            raise ValueError(f"Invalid --k_range: {s!r}. Must be positive and a<=b.")
        return ("range", a, b)

    k = int(s)
    if k <= 0:
        raise ValueError(f"Invalid --k_range: {s!r}. Must be a positive integer.")
    return ("single_k", k, k)


def read_input_lines(path: Path) -> List[Tuple[str, List[str]]]:
    """
    Reads the main mutations file:
      <pdb> <mut1,mut2,...>

    Ignores empty lines and lines starting with '#'.
    """
    rows: List[Tuple[str, List[str]]] = []
    with path.open("r", encoding="utf-8") as f:
        for i, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"Line {i} in {path} doesn't look like '<pdb> <mut1,mut2,...>': {raw!r}"
                )
            pdb = parts[0]
            mut_str = parts[1]
            muts = [m.strip() for m in mut_str.split(",") if m.strip()]
            if not muts:
                raise ValueError(f"Line {i} has no mutations after parsing: {raw!r}")
            rows.append((pdb, muts))
    if not rows:
        raise ValueError(f"No usable lines found in {path}")
    return rows


def enumerate_all_subsets(rows: Sequence[Tuple[str, List[str]]]) -> Dict[int, Set[MutCombo]]:
    """
    For each row (pdb, muts), generate all subsets of size 1..len(muts).
    Uniqueness is by (pdb, subset) where subset is order-insensitive.
    """
    by_k: Dict[int, Set[MutCombo]] = {}
    for pdb, muts in rows:
        uniq = sorted(set(muts))
        for k in range(1, len(uniq) + 1):
            by_k.setdefault(k, set())
            for combo in combinations(uniq, k):
                by_k[k].add(MutCombo(pdb=pdb, muts=tuple(sorted(combo))))
    return by_k


def parse_combo_line(raw: str) -> Optional[MutCombo]:
    """
    Parse a line from a completed file.
    Accepts:
      pdb<TAB>mut1,mut2
    or
      pdb mut1,mut2
    """
    line = raw.strip()
    if not line or line.startswith("#"):
        return None

    if "\t" in line:
        parts = line.split("\t")
    else:
        parts = line.split()

    if len(parts) < 2:
        return None

    pdb = parts[0].strip()
    mut_str = parts[1].strip()
    muts = tuple(sorted(m.strip() for m in mut_str.split(",") if m.strip()))
    if not pdb or not muts:
        return None
    return MutCombo(pdb=pdb, muts=muts)


def load_completed_files(paths: Sequence[Path]) -> Set[MutCombo]:
    completed: Set[MutCombo] = set()
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"--completed file not found: {p}")
        with p.open("r", encoding="utf-8") as f:
            for raw in f:
                mc = parse_combo_line(raw)
                if mc is not None:
                    completed.add(mc)
    return completed


def ks_to_process(mode: str, a: int, b: int, max_k: int) -> List[int]:
    if max_k <= 0:
        return []
    if mode == "all":
        return list(range(1, max_k + 1))
    start = max(1, a)
    end = min(max_k, b)
    if end < start:
        return []
    return list(range(start, end + 1))


def compute_progress_by_k(
    by_k: Dict[int, Set[MutCombo]],
    completed: Set[MutCombo],
    ks: Sequence[int],
) -> Dict[str, Dict[str, float]]:
    """
    Returns a dict keyed by str(k) with:
      total, completed, remaining, percent_completed
    Percent is 0 if total==0.
    """
    out: Dict[str, Dict[str, float]] = {}
    for k in ks:
        universe = by_k.get(k, set())
        total = len(universe)
        done = len(universe.intersection(completed)) if total else 0
        remaining = total - done
        pct = (done / total * 100.0) if total else 0.0
        out[str(k)] = {
            "total": float(total),
            "completed": float(done),
            "remaining": float(remaining),
            "percent_completed": float(pct),
        }
    return out


def select_combos(
    by_k: Dict[int, Set[MutCombo]],
    *,
    k_mode: str,
    k_a: int,
    k_b: int,
    n: int,
    random_sample: bool,
    seed: int,
    completed: Optional[Set[MutCombo]] = None,
) -> Dict[int, List[MutCombo]]:
    """
    Selection rules:
    - Choose which k to include via (k_mode, k_a, k_b). If k_mode=="all" => 1..max_k.
    - For each selected k:
        - Exclude combos in `completed` (if provided)
        - If n<=0 => return ALL remaining combos for that k (sorted)
        - Else => return up to n combos:
            - if random_sample: random sample without replacement (seeded), then sorted
            - else: first n from sorted list
    """
    if not by_k:
        return {}

    max_k = max(by_k.keys())
    ks = ks_to_process(k_mode, k_a, k_b, max_k)

    completed = completed or set()
    rng = random.Random(seed)

    out: Dict[int, List[MutCombo]] = {}
    for k in ks:
        all_list = sorted(by_k.get(k, set()))
        remaining = [c for c in all_list if c not in completed]

        if n <= 0:
            out[k] = remaining
            continue

        if not remaining:
            out[k] = []
            continue

        take = min(n, len(remaining))
        if random_sample:
            chosen = rng.sample(remaining, k=take)
            out[k] = sorted(chosen)
        else:
            out[k] = remaining[:take]

    return out


# =========================
# CLI / IO
# =========================

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_list(path: Path, combos: Sequence[MutCombo]) -> None:
    with path.open("w", encoding="utf-8") as f:
        for c in combos:
            f.write(c.to_line() + "\n")


def _tqdm(iterable: Iterable[int], **kwargs):
    if tqdm is None:
        return iterable
    return tqdm(iterable, **kwargs)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Enumerate mutation subsets (singles/doubles/...) and write selected outputs to a timestamped directory.",
        fromfile_prefix_chars="@",
    )
    ap.add_argument("--input_file", required=True, type=str, help="Path to mutations.txt-like input.")
    ap.add_argument(
        "--parent_out_dir",
        required=True,
        type=str,
        help="Parent output directory. A new timestamped subdir will be created inside.",
    )
    ap.add_argument("--run_name", required=True, type=str, help="Run name prefix for the timestamped output dir.")

    ap.add_argument(
        "--k_range",
        default="",
        type=str,
        help="Which mutation counts to output (e.g. 2-3). If omitted/empty => all k.",
    )

    ap.add_argument(
        "--n",
        default=0,
        type=int,
        help="How many to output per k (after excluding completed). 0 means 'all'.",
    )

    ap.add_argument(
        "--random",
        action="store_true",
        help="If set, randomly sample when --n > 0. Otherwise take first N sorted.",
    )
    ap.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed used only when --random is set (default: 0).",
    )

    ap.add_argument(
        "--completed",
        nargs="*",
        default=[],
        help="0+ txt files listing already-completed combos to exclude from selection.",
    )

    args = ap.parse_args()

    input_path = Path(args.input_file).expanduser().resolve()
    parent_out = Path(args.parent_out_dir).expanduser().resolve()

    k_mode, k_a, k_b = parse_k_range_arg(args.k_range)
    n = int(args.n)

    completed_paths = [Path(p).expanduser().resolve() for p in (args.completed or [])]
    completed_set = load_completed_files(completed_paths) if completed_paths else set()

    rows = read_input_lines(input_path)
    by_k = enumerate_all_subsets(rows)

    max_k = max(by_k.keys()) if by_k else 0
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = parent_out / f"{args.run_name}_{timestamp}"
    all_dir = out_dir / "all"
    sel_dir = out_dir / "selected"
    ensure_dir(all_dir)
    ensure_dir(sel_dir)

    all_ks = list(range(1, max_k + 1))
    ks_in_range = ks_to_process(k_mode, k_a, k_b, max_k)

    # NEW: progress summaries
    progress_all_by_k = compute_progress_by_k(by_k, completed_set, all_ks)
    progress_in_range_by_k = compute_progress_by_k(by_k, completed_set, ks_in_range)

    # summary
    summary = {
        "input_file": str(input_path),
        "parent_out_dir": str(parent_out),
        "run_name": args.run_name,
        "timestamp": timestamp,
        "max_k": max_k,
        "counts_all_by_k": {str(k): len(v) for k, v in sorted(by_k.items())},
        "k_range": args.k_range,
        "k_mode": k_mode,
        "n": n,
        "random": bool(args.random),
        "seed": int(args.seed),
        "completed_files": [str(p) for p in completed_paths],
        "completed_count_total": len(completed_set),
        "progress_all_by_k": progress_all_by_k,
        "progress_in_k_range_by_k": progress_in_range_by_k,
    }
    with (out_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Write ALL combos for every k
    for k in _tqdm(all_ks, desc="Writing ALL by k", unit="k"):
        write_list(all_dir / f"k{k}.txt", sorted(by_k.get(k, set())))

    # Selected with exclusion + N sampling
    selected = select_combos(
        by_k,
        k_mode=k_mode,
        k_a=k_a,
        k_b=k_b,
        n=n,
        random_sample=bool(args.random),
        seed=int(args.seed),
        completed=completed_set,
    )

    ks_selected = sorted(selected.keys())
    for k in _tqdm(ks_selected, desc="Writing SELECTED by k", unit="k"):
        write_list(sel_dir / f"k{k}.txt", selected[k])

    print(f"✅ Wrote outputs to: {out_dir}")
    print(f"   - All subsets:      {all_dir}")
    print(f"   - Selected subsets: {sel_dir}")
    print(f"   - Summary:          {out_dir / 'summary.json'}")


# =========================
# Unit tests
# =========================

class TestParseKRangeArg(unittest.TestCase):
    def test_empty_means_all(self):
        self.assertEqual(parse_k_range_arg(None)[0], "all")
        self.assertEqual(parse_k_range_arg("")[0], "all")
        self.assertEqual(parse_k_range_arg("   ")[0], "all")

    def test_range(self):
        self.assertEqual(parse_k_range_arg("2-3"), ("range", 2, 3))
        self.assertEqual(parse_k_range_arg("2-2"), ("range", 2, 2))

    def test_single_k(self):
        self.assertEqual(parse_k_range_arg("2"), ("single_k", 2, 2))

    def test_invalid(self):
        with self.assertRaises(ValueError):
            parse_k_range_arg("0")
        with self.assertRaises(ValueError):
            parse_k_range_arg("3-2")
        with self.assertRaises(ValueError):
            parse_k_range_arg("2-0")
        with self.assertRaises(ValueError):
            parse_k_range_arg("2-3-4")


class TestEnumerateAllSubsets(unittest.TestCase):
    def test_subset_counts_one_line_five_muts(self):
        rows = [("X.pdb", ["E19S", "S91K", "P149Y", "R214D", "G259R"])]
        by_k = enumerate_all_subsets(rows)
        self.assertEqual(len(by_k[1]), 5)
        self.assertEqual(len(by_k[2]), 10)
        self.assertEqual(len(by_k[3]), 10)
        self.assertEqual(len(by_k[4]), 5)
        self.assertEqual(len(by_k[5]), 1)

    def test_order_insensitive_and_dedup_within_line(self):
        rows = [("X.pdb", ["A1B", "C2D", "A1B", "C2D"])]
        by_k = enumerate_all_subsets(rows)
        self.assertEqual(len(by_k[1]), 2)
        self.assertEqual(len(by_k[2]), 1)

    def test_uniqueness_is_per_pdb(self):
        rows = [
            ("X.pdb", ["A1B", "C2D"]),
            ("Y.pdb", ["A1B", "C2D"]),
        ]
        by_k = enumerate_all_subsets(rows)
        self.assertEqual(len(by_k[2]), 2)


class TestSelection(unittest.TestCase):
    def setUp(self):
        self.rows = [("X.pdb", ["A", "B", "C", "D"])]
        self.by_k = enumerate_all_subsets(self.rows)

    def test_all_k_range_returns_all_ks(self):
        selected = select_combos(
            self.by_k,
            k_mode="all",
            k_a=1,
            k_b=1,
            n=0,
            random_sample=False,
            seed=0,
            completed=set(),
        )
        self.assertEqual(set(selected.keys()), {1, 2, 3, 4})
        self.assertEqual(len(selected[2]), 6)

    def test_range_only_doubles(self):
        selected = select_combos(
            self.by_k,
            k_mode="range",
            k_a=2,
            k_b=2,
            n=0,
            random_sample=False,
            seed=0,
            completed=set(),
        )
        self.assertEqual(set(selected.keys()), {2})
        self.assertEqual(len(selected[2]), 6)

    def test_n_limits_per_k(self):
        selected = select_combos(
            self.by_k,
            k_mode="all",
            k_a=1,
            k_b=1,
            n=2,
            random_sample=False,
            seed=0,
            completed=set(),
        )
        self.assertEqual(len(selected[1]), 2)
        self.assertEqual(len(selected[2]), 2)
        self.assertEqual(len(selected[3]), 2)
        self.assertEqual(len(selected[4]), 1)

    def test_random_is_deterministic_with_seed(self):
        s1 = select_combos(
            self.by_k,
            k_mode="all",
            k_a=1,
            k_b=1,
            n=3,
            random_sample=True,
            seed=123,
            completed=set(),
        )
        s2 = select_combos(
            self.by_k,
            k_mode="all",
            k_a=1,
            k_b=1,
            n=3,
            random_sample=True,
            seed=123,
            completed=set(),
        )
        self.assertEqual(s1, s2)

    def test_completed_exclusion_reduces_pool(self):
        completed = {MutCombo("X.pdb", ("A", "B"))}
        selected = select_combos(
            self.by_k,
            k_mode="range",
            k_a=2,
            k_b=2,
            n=0,
            random_sample=False,
            seed=0,
            completed=completed,
        )
        self.assertEqual(len(selected[2]), 5)
        self.assertNotIn(MutCombo("X.pdb", ("A", "B")), selected[2])

    def test_completed_exclusion_with_n_sampling(self):
        all_doubles = sorted(self.by_k[2])
        completed = set(all_doubles[:5])
        selected = select_combos(
            self.by_k,
            k_mode="range",
            k_a=2,
            k_b=2,
            n=50000,
            random_sample=True,
            seed=0,
            completed=completed,
        )
        self.assertEqual(len(selected[2]), 1)
        for c in selected[2]:
            self.assertNotIn(c, completed)


class TestCompletedParsing(unittest.TestCase):
    def test_parse_combo_line_tab(self):
        mc = parse_combo_line("X.pdb\tB,A\n")
        self.assertEqual(mc, MutCombo("X.pdb", ("A", "B")))

    def test_parse_combo_line_space(self):
        mc = parse_combo_line("X.pdb B,A\n")
        self.assertEqual(mc, MutCombo("X.pdb", ("A", "B")))

    def test_parse_combo_line_ignores_garbage(self):
        self.assertIsNone(parse_combo_line(""))
        self.assertIsNone(parse_combo_line("# comment"))
        self.assertIsNone(parse_combo_line("just_one_token"))


class TestProgressComputation(unittest.TestCase):
    def test_progress_counts(self):
        rows = [("X.pdb", ["A", "B", "C"])]  # k1=3, k2=3, k3=1
        by_k = enumerate_all_subsets(rows)

        completed = {
            MutCombo("X.pdb", ("A",)),
            MutCombo("X.pdb", ("B", "C")),
        }

        prog = compute_progress_by_k(by_k, completed, ks=[1, 2, 3])

        # k=1: total 3, completed 1, remaining 2
        self.assertEqual(prog["1"]["total"], 3.0)
        self.assertEqual(prog["1"]["completed"], 1.0)
        self.assertEqual(prog["1"]["remaining"], 2.0)

        # k=2: total 3, completed 1, remaining 2
        self.assertEqual(prog["2"]["total"], 3.0)
        self.assertEqual(prog["2"]["completed"], 1.0)
        self.assertEqual(prog["2"]["remaining"], 2.0)

        # k=3: total 1, completed 0, remaining 1
        self.assertEqual(prog["3"]["total"], 1.0)
        self.assertEqual(prog["3"]["completed"], 0.0)
        self.assertEqual(prog["3"]["remaining"], 1.0)

    def test_progress_percent_zero_total(self):
        # If ks includes a k not in by_k, total=0 and percent=0
        prog = compute_progress_by_k(by_k={}, completed=set(), ks=[5])
        self.assertEqual(prog["5"]["total"], 0.0)
        self.assertEqual(prog["5"]["percent_completed"], 0.0)


if __name__ == "__main__":
    main()