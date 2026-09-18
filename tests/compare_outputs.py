#!/usr/bin/env python3
"""
Compare pipeline outputs against a set of blessed expected outputs.

Usage
-----
Run pipeline first, then compare:
    python scripts/compare_outputs.py

Promote current test_output/ to the new expected baseline:
    python scripts/compare_outputs.py --bless

Override directories or seed:
    python scripts/compare_outputs.py --actual scripts/test_output --expected scripts/expected --seed 123
"""

import argparse
import filecmp
import os
import shutil
import sys

import pandas as pd

GREEN  = "\033[32m"
RED    = "\033[31m"
YELLOW = "\033[33m"
RESET  = "\033[0m"
BOLD   = "\033[1m"

SKIP_EXTENSIONS = {".log", ".trees"}


def _pass(name):
    print(f"  {GREEN}PASS{RESET}  {name}")


def _fail(name, detail):
    print(f"  {RED}FAIL{RESET}  {name}")
    for line in detail.splitlines():
        print(f"        {line}")


def _warn(name, detail):
    print(f"  {YELLOW}WARN{RESET}  {name}  —  {detail}")


# ---------------------------------------------------------------------------
# Per-type comparisons
# ---------------------------------------------------------------------------

MAX_DIFF_ROWS = 10


def compare_csv(name, exp_path, act_path):
    exp = pd.read_csv(exp_path)
    act = pd.read_csv(act_path)

    if exp.equals(act):
        _pass(name)
        return True

    lines = []

    # shape mismatch
    if exp.shape != act.shape:
        lines.append(f"shape: expected {exp.shape}, actual {act.shape}")

    # column mismatch
    if list(exp.columns) != list(act.columns):
        lines.append(f"columns: expected {list(exp.columns)}, actual {list(act.columns)}")
        _fail(name, "\n".join(lines))
        return False

    # same shape and columns: find which rows and cells differ
    if exp.shape == act.shape:
        row_mask = ~(exp == act).all(axis=1)
        n_diff_rows = row_mask.sum()
        lines.append(f"{n_diff_rows}/{len(exp)} row(s) differ")

        diff_rows = exp[row_mask].copy()
        diff_rows_act = act[row_mask].copy()
        for col in exp.columns:
            col_diff = exp.loc[row_mask, col] != act.loc[row_mask, col]
            if col_diff.any():
                lines.append(f"  column '{col}' differs in {col_diff.sum()} row(s):")
                sample = exp.loc[row_mask][col_diff][col].head(MAX_DIFF_ROWS)
                sample_act = act.loc[row_mask][col_diff][col].head(MAX_DIFF_ROWS)
                for idx, (e_val, a_val) in enumerate(zip(sample, sample_act)):
                    lines.append(f"    row {sample.index[idx]}: expected={e_val}  actual={a_val}")
                if col_diff.sum() > MAX_DIFF_ROWS:
                    lines.append(f"    ... and {col_diff.sum() - MAX_DIFF_ROWS} more")
    else:
        # different row counts: show row-level diff
        try:
            merged = exp.merge(act, indicator=True, how="outer")
            only_exp = merged[merged["_merge"] == "left_only"].drop(columns="_merge")
            only_act = merged[merged["_merge"] == "right_only"].drop(columns="_merge")
            if not only_exp.empty:
                lines.append(f"rows only in expected (first {MAX_DIFF_ROWS}):")
                lines.append(only_exp.head(MAX_DIFF_ROWS).to_string(index=False))
            if not only_act.empty:
                lines.append(f"rows only in actual (first {MAX_DIFF_ROWS}):")
                lines.append(only_act.head(MAX_DIFF_ROWS).to_string(index=False))
        except Exception as e:
            lines.append(f"could not diff rows: {e}")

    _fail(name, "\n".join(lines))
    return False


def compare_vcf(name, exp_path, act_path):
    """Compare VCF files ignoring ## meta-header lines (may contain tool versions)."""
    def data_lines(path):
        with open(path) as f:
            return [l for l in f if not l.startswith("##")]

    exp_lines = data_lines(exp_path)
    act_lines = data_lines(act_path)

    if exp_lines == act_lines:
        _pass(name)
        return True

    exp_set = set(exp_lines)
    act_set = set(act_lines)
    only_exp = exp_set - act_set
    only_act = act_set - exp_set

    lines = [f"{len(only_exp)} line(s) only in expected, {len(only_act)} line(s) only in actual"]
    if only_exp:
        lines.append("expected-only (first 5):")
        for l in sorted(only_exp)[:5]:
            lines.append("  " + l.rstrip())
    if only_act:
        lines.append("actual-only (first 5):")
        for l in sorted(only_act)[:5]:
            lines.append("  " + l.rstrip())

    _fail(name, "\n".join(lines))
    return False



# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

def compare_file(name, exp_path, act_path):
    ext = os.path.splitext(name)[1].lower()
    if ext == ".csv":
        return compare_csv(name, exp_path, act_path)
    elif ext == ".vcf":
        return compare_vcf(name, exp_path, act_path)
    else:
        if filecmp.cmp(exp_path, act_path, shallow=False):
            _pass(name)
            return True
        _fail(name, "files differ (binary comparison)")
        return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compare pipeline outputs against expected outputs")
    parser.add_argument("--expected", default="scripts/expected",
                        help="Directory of blessed expected outputs (default: scripts/expected)")
    parser.add_argument("--actual",   default="scripts/test_output",
                        help="Directory of outputs from the latest run (default: scripts/test_output)")
    parser.add_argument("--bless", action="store_true",
                        help="Copy actual/ -> expected/ to set a new baseline (skips log files)")
    args = parser.parse_args()

    if args.bless:
        os.makedirs(args.expected, exist_ok=True)
        blessed = []
        skipped = []
        for fname in sorted(os.listdir(args.actual)):
            if os.path.splitext(fname)[1].lower() in SKIP_EXTENSIONS:
                skipped.append(fname)
                continue
            src = os.path.join(args.actual, fname)
            dst = os.path.join(args.expected, fname)
            shutil.copy2(src, dst)
            blessed.append(fname)
        print(f"\nBlessed {len(blessed)} file(s) to {args.expected}/:")
        for f in blessed:
            print(f"  {f}")
        if skipped:
            print(f"Skipped {len(skipped)} log file(s): {', '.join(skipped)}")
        print(f"\nCommit {args.expected}/ to lock in the new baseline.")
        return

    if not os.path.isdir(args.expected):
        print(f"{RED}Expected directory not found: {args.expected}{RESET}")
        print("Run with --bless after a clean pipeline run to create it.")
        sys.exit(1)

    expected_files = {
        f for f in os.listdir(args.expected)
        if os.path.isfile(os.path.join(args.expected, f))
        and os.path.splitext(f)[1].lower() not in SKIP_EXTENSIONS
    }

    print(f"\n{BOLD}Comparing {args.actual}/ against {args.expected}/{RESET}\n")

    results = []

    for fname in sorted(expected_files):
        exp_path = os.path.join(args.expected, fname)
        act_path = os.path.join(args.actual, fname)
        if not os.path.exists(act_path):
            _warn(fname, f"missing from {args.actual}/")
            results.append(False)
            continue
        results.append(compare_file(fname, exp_path, act_path))

    # report files in actual that aren't in expected (new outputs, not yet blessed)
    actual_files = {
        f for f in os.listdir(args.actual)
        if os.path.isfile(os.path.join(args.actual, f))
        and os.path.splitext(f)[1].lower() not in SKIP_EXTENSIONS
    }
    new_files = actual_files - expected_files
    if new_files:
        print()
        for f in sorted(new_files):
            _warn(f, "new output file — not yet in expected/. Answer 'y' at the bless prompt below to add it.")

    n_pass = sum(results)
    n_fail = len(results) - n_pass
    print()
    if n_fail == 0:
        print(f"{GREEN}{BOLD}All {n_pass} check(s) passed.{RESET}")
        sys.exit(0)
    else:
        print(f"{RED}{BOLD}{n_fail}/{len(results)} check(s) failed.{RESET}")
        sys.exit(1)


if __name__ == "__main__":
    main()
