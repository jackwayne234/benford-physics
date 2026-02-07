#!/usr/bin/env python3
"""
run_test.py -- Sequential test runner for the Benford's Law testing lab.

Usage:
    python3 scripts/run_test.py <test_id>          Run one specific test
    python3 scripts/run_test.py --next              Run the next queued test
    python3 scripts/run_test.py --all               Run all queued tests sequentially
    python3 scripts/run_test.py --retry-errors      Retry all tests with 'error' status
    python3 scripts/run_test.py --status            Show queue status summary
"""

import argparse
import fcntl
import json
import os
import sys
import traceback
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

sys.path.insert(0, SCRIPT_DIR)

from benford_core import run_full_analysis
from data_fetchers import FETCHER_REGISTRY

QUEUE_PATH = os.path.join(PROJECT_ROOT, "data", "test_queue.json")
SOURCES_DIR = os.path.join(PROJECT_ROOT, "data", "sources")
INDIVIDUAL_DIR = os.path.join(PROJECT_ROOT, "results", "individual")
SUMMARY_PATH = os.path.join(PROJECT_ROOT, "results", "summary.json")

# ---------------------------------------------------------------------------
# File I/O helpers (with locking)
# ---------------------------------------------------------------------------


def load_json_locked(path):
    """Read a JSON file under a shared (read) lock."""
    with open(path, "r") as f:
        fcntl.flock(f, fcntl.LOCK_SH)
        try:
            data = json.load(f)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)
    return data


def write_json_locked(path, data):
    """Write a JSON file under an exclusive lock."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            json.dump(data, f, indent=2)
            f.write("\n")
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def load_queue():
    """Load the test queue from disk. Wraps plain arrays into {'tests': [...]}."""
    data = load_json_locked(QUEUE_PATH)
    if isinstance(data, list):
        return {"tests": data}
    return data


def save_queue(queue):
    """Save the test queue back to disk. Unwraps to plain array."""
    data = queue.get("tests", queue) if isinstance(queue, dict) else queue
    write_json_locked(QUEUE_PATH, data)


# ---------------------------------------------------------------------------
# Queue manipulation
# ---------------------------------------------------------------------------


def find_test_entry(queue, test_id):
    """Return (index, entry) for a test_id, or (None, None)."""
    for i, entry in enumerate(queue.get("tests", [])):
        if entry.get("test_id") == test_id:
            return i, entry
    return None, None


def set_test_status(test_id, status, error_message=None):
    """Atomically update a test's status in the queue file."""
    queue = load_queue()
    idx, entry = find_test_entry(queue, test_id)
    if idx is not None:
        queue["tests"][idx]["status"] = status
        if error_message:
            queue["tests"][idx]["error_message"] = error_message
        elif "error_message" in queue["tests"][idx] and status != "error":
            del queue["tests"][idx]["error_message"]
        save_queue(queue)


def get_next_queued(queue, status_filter="queued"):
    """Return the next test entry with the given status, respecting priority."""
    candidates = [
        e for e in queue.get("tests", []) if e.get("status") == status_filter
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda e: (e.get("priority", 999), 0))
    return candidates[0]


# ---------------------------------------------------------------------------
# Core test execution
# ---------------------------------------------------------------------------


def execute_test(test_id):
    """Run a single Benford's Law test end-to-end. Returns True on success."""

    # Step 1-2: Load queue and find entry
    queue = load_queue()
    idx, test_entry = find_test_entry(queue, test_id)
    if idx is None:
        print(f"  ERROR: test_id '{test_id}' not found in queue.")
        return False

    # Step 3: Mark in_progress
    set_test_status(test_id, "in_progress")

    try:
        # Step 4: Fetch data
        if test_id not in FETCHER_REGISTRY:
            raise KeyError(f"No fetcher registered for '{test_id}'")
        fetch_result = FETCHER_REGISTRY[test_id]()
        # Fetcher returns a dict with 'values' key (list of floats)
        if isinstance(fetch_result, dict):
            values = fetch_result.get("values", [])
        else:
            values = list(fetch_result)

        # Step 5: Validate minimum data points
        if len(values) < 50:
            set_test_status(test_id, "insufficient_data")
            print(f"  SKIPPED: Only {len(values)} data points (need >= 50).")
            return False

        # Step 6: Cache raw data
        os.makedirs(SOURCES_DIR, exist_ok=True)
        cache_obj = {
            "test_id": test_id,
            "fetch_timestamp": datetime.now(timezone.utc).isoformat(),
            "raw_count": fetch_result.get("raw_count", len(values)) if isinstance(fetch_result, dict) else len(values),
            "valid_count": len(values),
            "values": values,
            "notes": fetch_result.get("notes", "") if isinstance(fetch_result, dict) else "",
        }
        cache_path = os.path.join(SOURCES_DIR, f"{test_id}.json")
        write_json_locked(cache_path, cache_obj)

        # Step 7: Run analysis
        analysis = run_full_analysis(values)

        # Step 8: Build result dict
        fd = analysis.get("first_digit", {})
        result = {
            "test_id": test_id,
            "display_name": test_entry.get("display_name", test_id),
            "category": test_entry.get("category", ""),
            "description": test_entry.get("description", ""),
            "status": "complete",
            "data_points": analysis.get("n", 0),
            "source": test_entry.get("data_source", ""),
            "digit_distribution": fd.get("observed_distribution", {}),
            "expected_distribution": fd.get("expected_distribution", {}),
            "per_digit_deviation": fd.get("per_digit_deviation", {}),
            "chi_squared": fd.get("chi_squared", {}),
            "mad": fd.get("mad", {}).get("mad", 0.0),
            "mad_classification": fd.get("mad", {}).get("classification", ""),
            "ks_test": fd.get("ks_test", {}),
            "delta_b": fd.get("euclidean_deviation", 0.0),
            "verdict": analysis.get("verdict", "MARGINAL"),
            "confidence": "high" if analysis.get("n", 0) >= 100 else "low",
            "notes": analysis.get("notes", ""),
            "interesting_findings": "",
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }

        # Step 9: Write individual result
        os.makedirs(INDIVIDUAL_DIR, exist_ok=True)
        result_path = os.path.join(INDIVIDUAL_DIR, f"{test_id}.json")
        write_json_locked(result_path, result)

        # Step 10: Update summary.json
        update_summary(result)

        # Step 11: Mark complete
        set_test_status(test_id, "complete")

        # Step 12: Print summary
        delta_b = result["delta_b"]
        verdict = result["verdict"]
        dp = result["data_points"]
        conf = result["confidence"]
        symbol = "+" if verdict == "CONFORMS" else ("-" if verdict == "DEVIATES" else "~")
        print(f"  {verdict} (d_B={delta_b:.3f}, n={dp}, confidence={conf}) [{symbol}]")
        return True

    except Exception as exc:
        set_test_status(test_id, "error", error_message=str(exc))
        print(f"  ERROR: {exc}")
        traceback.print_exc()
        return False


# ---------------------------------------------------------------------------
# Summary update
# ---------------------------------------------------------------------------


def update_summary(result):
    """Load summary.json, update/append the entry, recalculate stats, write back."""

    # Load or initialize
    if os.path.exists(SUMMARY_PATH):
        summary = load_json_locked(SUMMARY_PATH)
    else:
        summary = {"tests": [], "stats": {}, "last_updated": ""}

    if "tests" not in summary:
        summary["tests"] = []
    if "stats" not in summary:
        summary["stats"] = {}

    # Build summary entry
    chi_p = result.get("chi_squared", {}).get("p_value", None)
    ks_p = result.get("ks_test", {}).get("p_value", None)
    entry = {
        "test_id": result["test_id"],
        "display_name": result["display_name"],
        "category": result["category"],
        "data_points": result["data_points"],
        "verdict": result["verdict"],
        "delta_b": result["delta_b"],
        "mad": result["mad"],
        "chi_squared_p": chi_p,
        "ks_p": ks_p,
        "confidence": result["confidence"],
        "timestamp": result["timestamp"],
    }

    # Find or create entry
    found = False
    for i, existing in enumerate(summary["tests"]):
        if existing.get("test_id") == result["test_id"]:
            summary["tests"][i] = entry
            found = True
            break
    if not found:
        summary["tests"].append(entry)

    # Recalculate stats
    completed = [t for t in summary["tests"] if t.get("verdict")]
    summary["stats"] = {
        "total_completed": len(completed),
        "conforms": sum(1 for t in completed if t["verdict"] == "CONFORMS"),
        "marginal": sum(1 for t in completed if t["verdict"] == "MARGINAL"),
        "deviates": sum(1 for t in completed if t["verdict"] == "DEVIATES"),
        "interesting": sum(1 for t in completed if t["verdict"] == "INTERESTING"),
    }
    summary["last_updated"] = datetime.now(timezone.utc).isoformat()

    write_json_locked(SUMMARY_PATH, summary)


# ---------------------------------------------------------------------------
# --status command
# ---------------------------------------------------------------------------


def show_status():
    """Print a formatted queue status summary."""
    queue = load_queue()
    tests = queue.get("tests", [])
    total = len(tests)

    status_counts = {}
    category_data = {}

    for t in tests:
        st = t.get("status", "unknown")
        status_counts[st] = status_counts.get(st, 0) + 1

        cat = t.get("category", "Uncategorized")
        if cat not in category_data:
            category_data[cat] = {"total": 0, "complete": 0}
        category_data[cat]["total"] += 1
        if st == "complete":
            category_data[cat]["complete"] += 1

    print()
    print("Benford's Law Testing Lab -- Queue Status")
    print("=" * 42)
    print(f"Total tests: {total}")
    print(f"Queued:      {status_counts.get('queued', 0)}")
    print(f"In Progress: {status_counts.get('in_progress', 0)}")
    print(f"Complete:    {status_counts.get('complete', 0)}")
    print(f"Errors:      {status_counts.get('error', 0)}")
    insufficient = status_counts.get("insufficient_data", 0)
    if insufficient:
        print(f"Insufficient: {insufficient}")
    print()
    print("By Category:")

    max_cat_len = max((len(c) for c in category_data), default=10)
    for cat in sorted(category_data):
        info = category_data[cat]
        label = f"{cat}:".ljust(max_cat_len + 2)
        print(f"  {label} {info['complete']}/{info['total']} complete")
    print()


# ---------------------------------------------------------------------------
# --all mode
# ---------------------------------------------------------------------------


def run_all_queued():
    """Run all queued tests in priority order."""
    queue = load_queue()
    candidates = [
        t for t in queue.get("tests", []) if t.get("status") == "queued"
    ]
    if not candidates:
        print("No queued tests remaining.")
        return

    # Sort by priority then queue order
    candidates.sort(key=lambda t: (t.get("priority", 999), 0))
    total_all = len(queue.get("tests", []))
    num_queued = len(candidates)

    print(f"\nRunning {num_queued} queued tests (out of {total_all} total)...\n")

    for i, entry in enumerate(candidates, 1):
        test_id = entry["test_id"]
        display = entry.get("display_name", test_id)
        print(f"[{i}/{num_queued}] Running: {test_id} ({display}) ...", end=" ", flush=True)
        success = execute_test(test_id)
        if not success:
            # Message already printed inside execute_test
            pass
    print(f"\nDone. Ran {num_queued} tests.\n")


# ---------------------------------------------------------------------------
# --retry-errors mode
# ---------------------------------------------------------------------------


def retry_errors():
    """Re-run all tests that previously failed with 'error' status."""
    queue = load_queue()
    error_tests = [
        t for t in queue.get("tests", []) if t.get("status") == "error"
    ]
    if not error_tests:
        print("No failed tests to retry.")
        return

    num_errors = len(error_tests)
    print(f"\nRetrying {num_errors} failed test(s)...\n")

    for i, entry in enumerate(error_tests, 1):
        test_id = entry["test_id"]
        display = entry.get("display_name", test_id)
        # Reset status to queued before retrying
        set_test_status(test_id, "queued")
        print(f"[{i}/{num_errors}] Retrying: {test_id} ({display}) ...", end=" ", flush=True)
        execute_test(test_id)

    print(f"\nDone. Retried {num_errors} tests.\n")


# ---------------------------------------------------------------------------
# --next mode
# ---------------------------------------------------------------------------


def run_next():
    """Run the single next queued test."""
    queue = load_queue()
    entry = get_next_queued(queue)
    if entry is None:
        print("No queued tests remaining.")
        return

    test_id = entry["test_id"]
    display = entry.get("display_name", test_id)
    print(f"\nRunning next queued test: {test_id} ({display}) ...", end=" ", flush=True)
    execute_test(test_id)
    print()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Benford's Law Testing Lab -- Sequential Test Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/run_test.py bose_einstein_numerical\n"
            "  python3 scripts/run_test.py --next\n"
            "  python3 scripts/run_test.py --all\n"
            "  python3 scripts/run_test.py --retry-errors\n"
            "  python3 scripts/run_test.py --status\n"
        ),
    )
    parser.add_argument(
        "test_id",
        nargs="?",
        default=None,
        help="Specific test_id to run",
    )
    parser.add_argument(
        "--next",
        action="store_true",
        dest="run_next",
        help="Run the next queued test",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        dest="run_all",
        help="Run all queued tests sequentially",
    )
    parser.add_argument(
        "--retry-errors",
        action="store_true",
        dest="retry_errors",
        help="Retry all tests with 'error' status",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Show queue status summary",
    )

    args = parser.parse_args()

    # Determine which mode to run
    mode_count = sum([
        args.test_id is not None,
        args.run_next,
        args.run_all,
        args.retry_errors,
        args.status,
    ])

    if mode_count == 0:
        parser.print_help()
        sys.exit(1)

    if mode_count > 1:
        print("ERROR: Please specify only one mode at a time.")
        sys.exit(1)

    # Validate that the queue file exists (required for all modes)
    if not os.path.exists(QUEUE_PATH):
        print(f"ERROR: Queue file not found: {QUEUE_PATH}")
        sys.exit(1)

    # Dispatch
    if args.status:
        show_status()
    elif args.run_all:
        run_all_queued()
    elif args.run_next:
        run_next()
    elif args.retry_errors:
        retry_errors()
    elif args.test_id:
        test_id = args.test_id
        queue = load_queue()
        _, entry = find_test_entry(queue, test_id)
        if entry is None:
            print(f"ERROR: test_id '{test_id}' not found in queue.")
            sys.exit(1)
        display = entry.get("display_name", test_id)
        print(f"\nRunning: {test_id} ({display}) ...", end=" ", flush=True)
        success = execute_test(test_id)
        if not success:
            sys.exit(1)
        print()


if __name__ == "__main__":
    main()
