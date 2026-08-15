"""Fail-fast verification for the high-school RNALocate analysis package."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PAPER = ROOT

EXPECTED_MODELS = {
    "prevalence",
    "context_only",
    "support_only",
    "sequence_only",
    "context_sequence",
    "context_support",
    "support_sequence",
    "full",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def close(actual: float, expected: float, tolerance: float = 1e-9) -> None:
    require(abs(actual - expected) <= tolerance, f"Expected {expected}, found {actual}")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def verify_metrics() -> None:
    metrics = pd.read_csv(RESULTS / "record_availability_metrics.tsv", sep="\t")
    require(set(metrics["model"]) == EXPECTED_MODELS, "Factorial model matrix is incomplete")
    require(len(metrics) == 8, "Expected exactly eight factorial models")
    lookup = metrics.set_index("model")
    close(float(lookup.loc["context_sequence", "average_precision"]), 0.8823216328735735)
    close(float(lookup.loc["support_sequence", "average_precision"]), 0.9071274953008045)
    close(float(lookup.loc["context_support", "average_precision"]), 0.9571531912808816)
    close(float(lookup.loc["full", "average_precision"]), 0.9580277680595601)
    require((metrics["ap_ci95_low"] <= metrics["average_precision"]).all(), "AP below lower interval")
    require((metrics["average_precision"] <= metrics["ap_ci95_high"]).all(), "AP above upper interval")

    effects = pd.read_csv(RESULTS / "record_availability_incremental_effects.tsv", sep="\t")
    require(len(effects) == 4, "Expected four prespecified incremental effects")
    for row in effects.itertuples(index=False):
        expected = float(lookup.loc[row.left_model, "average_precision"] - lookup.loc[row.right_model, "average_precision"])
        close(float(row.ap_difference), expected)
        require(row.ap_difference_ci95_low <= row.ap_difference <= row.ap_difference_ci95_high, f"Effect interval misses {row.contrast}")


def verify_secondary_analyses() -> None:
    gradient = pd.read_csv(RESULTS / "support_gradient.tsv", sep="\t")
    require(gradient["other_context_support"].tolist() == list(range(16)), "Support bins must be 0 through 15")
    require(int(gradient["test_pairs"].sum()) == 45152, "Support gradient does not cover the full test grid")
    zero = gradient.loc[gradient["other_context_support"].eq(0)].iloc[0]
    close(float(zero["observed_fraction"]), 1.0)
    require(bool(zero["selection_boundary"]), "Support zero must be marked as a selection boundary")

    permutation = pd.read_csv(RESULTS / "shuffled_support_control.tsv", sep="\t")
    require(len(permutation) == 200, "Expected 200 support permutations")
    require(permutation["seed"].nunique() == 200, "Permutation seeds must be unique")
    require(float(permutation["average_precision"].max()) < float(permutation["real_context_support_ap"].min()), "Null control overlaps the real AP")

    context = pd.read_csv(RESULTS / "record_availability_by_context.tsv", sep="\t")
    require(context["context"].nunique() == 16, "Expected 16 contexts")
    require(len(context) == 64, "Expected four model rows for each context")
    require(np.isfinite(context[["average_precision", "brier", "ap_lift"]].to_numpy()).all(), "Non-finite context metric")

    calibration = pd.read_csv(RESULTS / "calibration_bins.tsv", sep="\t")
    for model, frame in calibration.groupby("model"):
        require(int(frame["test_pairs"].sum()) == 45152, f"Calibration bins do not cover test rows for {model}")
        require(((frame["observed_fraction"] >= 0) & (frame["observed_fraction"] <= 1)).all(), "Invalid calibration fraction")

    sensitivity = pd.read_csv(RESULTS / "sequence_view_sensitivity.tsv", sep="\t")
    require(len(sensitivity) == 6, "Expected six sequence-view sensitivity rows")
    require(set(sensitivity["feature_view"]) == {"composition_length", "three_mer_only", "complete_70"}, "Missing sequence view")

    audit = json.loads((RESULTS / "audit_report.json").read_text(encoding="utf-8"))
    require(audit["status"] == "pass", "Analysis audit did not pass")
    require(audit["support_zero_is_all_observed"] is True, "Support-zero boundary audit failed")

    support_zero = pd.read_csv(RESULTS / "support_zero_exclusion_metrics.tsv", sep="\t").set_index("model")
    require(set(support_zero.index) == {"context_only", "support_only", "context_support", "full"}, "Support-zero model set is incomplete")
    close(float(support_zero.loc["context_only", "average_precision"]), 0.8014217430297013)
    close(float(support_zero.loc["support_only", "average_precision"]), 0.9005935457628997)
    close(float(support_zero.loc["context_support", "average_precision"]), 0.9573395950016954)
    close(float(support_zero.loc["full", "average_precision"]), 0.9581977967677284)
    require((support_zero["test_pairs"] == 45003).all(), "Support-zero sensitivity test size changed")
    support_zero_audit = json.loads((RESULTS / "support_zero_exclusion_audit.json").read_text(encoding="utf-8"))
    require(support_zero_audit["status"] == "pass", "Support-zero sensitivity audit did not pass")
    require(support_zero_audit["excluded_pairs"] == 992, "Unexpected support-zero exclusion count")
    require(support_zero_audit["retained_pairs"] + support_zero_audit["excluded_pairs"] == support_zero_audit["full_pairs"], "Support-zero pair accounting failed")


def verify_figures() -> None:
    manifest = json.loads((PAPER / "figures/figure_manifest.json").read_text(encoding="utf-8"))
    require(len(manifest["files"]) == 12, "Expected four figures in three formats")
    require(len(manifest["source_groups"]) == 4, "Expected four source-data groups")
    for group in manifest["source_groups"].values():
        for relative in group:
            require((ROOT / relative).exists(), f"Missing source table {relative}")
    for entry in manifest["files"]:
        path = ROOT / entry["file"]
        require(path.exists() and path.stat().st_size > 1000, f"Missing or empty figure {path}")
        require(digest(path) == entry["sha256"], f"Hash mismatch for {path}")
        if entry["format"] == "png":
            with Image.open(path) as image:
                require(image.mode == "RGB", f"PNG must have an opaque RGB background: {path}")
                require(image.width >= 2100 and image.height >= 2000, f"PNG dimensions are too small: {path}")
                dpi = image.info.get("dpi", (0, 0))
                require(min(dpi) >= 299, f"PNG DPI metadata below 299: {path}")
    captions = (PAPER / "figure_captions.md").read_text(encoding="utf-8")
    for number in range(1, 5):
        require(f"## Figure {number}" in captions, f"Missing Figure {number} caption")
    require(captions.count("**Alt text.**") == 4, "Every figure requires alt text")


def verify_all() -> None:
    verify_metrics()
    verify_secondary_analyses()
    verify_figures()
    for relative in ["paper/preprint.pdf", "paper/preprint.docx", "paper/supporting_information.pdf", "paper/manuscript.md"]:
        require((ROOT / relative).exists() and (ROOT / relative).stat().st_size > 1000, f"Missing paper artifact: {relative}")


def main() -> None:
    verify_all()
    print("PASS: high-school context-shortcut analysis package verified")


if __name__ == "__main__":
    main()
