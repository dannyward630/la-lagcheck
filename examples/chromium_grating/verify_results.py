"""Compare completed grating runs with the packaged numerical expectations.

Exit 0: matched; exit 2: mismatch; exit 1: unreadable input or tool error.
This checks the saved outputs, not geometric accuracy or an importer round trip.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(actual, expected, label, mismatches):
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            mismatches.append(label + ": expected an object")
            return
        for key, value in expected.items():
            if key not in actual:
                mismatches.append(label + "." + key + ": missing")
            else:
                compare(actual[key], value, label + "." + key, mismatches)
    elif isinstance(expected, bool) or expected is None:
        if actual is not expected:
            mismatches.append(label + ": differs")
    elif isinstance(expected, (int, float)):
        if (isinstance(actual, bool) or not isinstance(actual, (int, float))
                or not math.isfinite(actual)
                or not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-6)):
            mismatches.append(label + ": numeric value differs")
    elif actual != expected:
        mismatches.append(label + ": differs")


def path_value(data, path):
    for key in path.split("."):
        data = data[key]
    return data


def verify(automatic, declared):
    expected = read_json(HERE / "expected_results.json")
    mismatches = []
    six_checks = {}
    for name, folder in [("automatic", automatic), ("declared", declared)]:
        report = read_json(folder / "report.json")
        compare(report, expected["runs"][name]["report_subset"], name, mismatches)
        compare(report.get("input_sha256"), expected["input_sha256"], name + ".input_sha256", mismatches)
        six_checks[name] = {}
        for group, fields in expected["six_check_report_fields"].items():
            try:
                passed = all(path_value(report, field) is True for field in fields)
            except (KeyError, TypeError):
                passed = False
            six_checks[name][group] = passed
            if not passed:
                mismatches.append(name + ".six_checks." + group + ": not met")
        if (folder / "INCOMPLETE.json").exists():
            mismatches.append(name + ": INCOMPLETE.json exists")
        if name == "declared":
            post = read_json(folder / "post_export_report.json")
            compare(post, expected["runs"][name]["post_export_report_subset"], "declared.post_export_report", mismatches)
            manifest = read_json(folder / "correction_manifest.json")
            compare(manifest, expected["runs"][name]["manifest_subset"], "declared.manifest", mismatches)
            candidate = sha256(folder / "candidate_laser_log.csv")
            corrected = sha256(folder / "corrected_laser_log.csv")
            recorded = report.get("post_export_verification", {}).get("diagnostic_input_sha256")
            if not (candidate == corrected == recorded == manifest.get("output_sha256")):
                mismatches.append("declared: candidate, final and diagnosed bytes differ")
    for label, relative in expected["input_files"].items():
        if sha256(HERE / relative) != expected["input_sha256"][label]:
            mismatches.append(label + ": packaged raw input bytes differ")
    return {"status": "matched" if not mismatches else "mismatch",
            "mismatches": mismatches, "six_checks": six_checks,
            "scope": "Saved numerical results and raw input bytes; no geometric or receiver validation"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--automatic", type=Path, required=True)
    parser.add_argument("--declared", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = verify(args.automatic, args.declared)
    except Exception as exc:
        print(json.dumps({"status": "tool_or_input_error", "error": type(exc).__name__ + ": " + str(exc)}))
        return 1
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "matched" else 2


if __name__ == "__main__":
    raise SystemExit(main())
