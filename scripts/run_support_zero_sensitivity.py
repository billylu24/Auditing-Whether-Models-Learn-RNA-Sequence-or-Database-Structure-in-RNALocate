"""Run the support-zero exclusion sensitivity analysis.

This post hoc sensitivity analysis was developed in response to external feedback. It removes every
gene-context pair with other-context support equal to zero from training,
validation, and test sets, then refits the four requested models.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_analysis import (
    WeightedAveragePrecision,
    build_grid,
    design_matrix,
    fit_model,
    safe_metrics,
)


MODELS: dict[str, tuple[str, ...]] = {
    "context_only": ("context",),
    "support_only": ("support",),
    "context_support": ("context", "support"),
    "full": ("context", "support", "sequence"),
}
DISPLAY_NAMES = {
    "context_only": "Context",
    "support_only": "Other-context support",
    "context_support": "Context + support",
    "full": "Context + support + sequence",
}
CONTRASTS = {
    "support_after_context": ("context_support", "context_only"),
    "sequence_after_context_support": ("full", "context_support"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=ROOT / "inputs/eligible_genes.tsv.gz")
    parser.add_argument("--reference", type=Path, default=ROOT / "inputs/observed_gene_context_pairs.tsv.gz")
    parser.add_argument("--split", type=Path, default=ROOT / "inputs/gene_component_split.tsv")
    parser.add_argument("--features", type=Path, default=ROOT / "inputs/sequence_features_70.npy")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260815)
    return parser.parse_args()


def bootstrap_models(predictions: pd.DataFrame, replicates: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    components, component_index = np.unique(predictions["cluster_representative"].to_numpy(), return_inverse=True)
    y = predictions["observed"].to_numpy(dtype=np.uint8)
    calculators = {
        model: WeightedAveragePrecision(y, predictions[f"prob_{model}"].to_numpy())
        for model in MODELS
    }
    values = np.empty((replicates, len(MODELS)), dtype=np.float64)
    rng = np.random.default_rng(seed)
    for replicate in range(replicates):
        sampled = rng.integers(0, len(components), size=len(components))
        component_weights = np.bincount(sampled, minlength=len(components))
        row_weights = component_weights[component_index]
        for column, calculator in enumerate(calculators.values()):
            values[replicate, column] = calculator.score(row_weights)

    bootstrap = pd.DataFrame(values, columns=list(MODELS))
    intervals = pd.DataFrame(
        [
            {
                "model": model,
                "ap_ci95_low": float(np.nanquantile(bootstrap[model], 0.025)),
                "ap_ci95_high": float(np.nanquantile(bootstrap[model], 0.975)),
                "bootstrap_replicates": replicates,
                "bootstrap_unit": "all-sequence graph component",
            }
            for model in MODELS
        ]
    )
    effects = []
    for contrast, (left, right) in CONTRASTS.items():
        difference = bootstrap[left] - bootstrap[right]
        effects.append(
            {
                "contrast": contrast,
                "left_model": left,
                "right_model": right,
                "ap_difference_ci95_low": float(np.nanquantile(difference, 0.025)),
                "ap_difference_ci95_high": float(np.nanquantile(difference, 0.975)),
                "bootstrap_replicates": replicates,
                "bootstrap_unit": "all-sequence graph component",
            }
        )
    return bootstrap, intervals, pd.DataFrame(effects)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    full_grid, cached = build_grid(args)
    reduced = full_grid.loc[full_grid["other_context_support"].gt(0)].copy().reset_index(drop=True)

    train = reduced["split"].eq("train").to_numpy()
    validation = reduced["split"].eq("validation").to_numpy()
    test = reduced["split"].eq("test").to_numpy()
    y = reduced["observed"].to_numpy(dtype=np.uint8)
    gene = reduced["gene_index"].to_numpy()
    context = reduced["context_index"].to_numpy()
    support = reduced["other_context_support"].to_numpy()
    context_count = int(context.max()) + 1

    predictions = reduced.loc[test, ["rna_symbol", "context", "cluster_representative", "observed", "other_context_support"]].reset_index(drop=True)
    metric_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    for model_name, parts in MODELS.items():
        x = design_matrix(parts, gene, context, support, cached, context_count)
        fitted = fit_model(x, y, train, validation, test, model_name, args.seed)
        predictions[f"prob_{model_name}"] = fitted.probability
        selection_rows.extend(fitted.selection_rows)
        metric_rows.append(
            {
                "analysis_set": "other_context_support_gt_0",
                "model": model_name,
                "display_name": DISPLAY_NAMES[model_name],
                "selected_C": fitted.selected_c,
                "validation_average_precision": fitted.validation_ap,
                **safe_metrics(y[test], fitted.probability),
            }
        )
        del x
        gc.collect()

    bootstrap, intervals, effects = bootstrap_models(predictions, args.bootstrap, args.seed + 1000)
    metrics = pd.DataFrame(metric_rows).merge(intervals, on="model", how="left", validate="one_to_one")
    lookup = metrics.set_index("model")["average_precision"]
    effects["ap_difference"] = [float(lookup.loc[left] - lookup.loc[right]) for left, right in zip(effects["left_model"], effects["right_model"])]

    split_counts = full_grid.groupby("split", observed=True).size().rename("full_pairs").to_frame()
    split_counts["retained_pairs"] = reduced.groupby("split", observed=True).size()
    split_counts["excluded_support_zero_pairs"] = split_counts["full_pairs"] - split_counts["retained_pairs"]
    split_counts = split_counts.reset_index()

    test_y = y[test]
    metrics.insert(1, "test_pairs", int(test.sum()))
    metrics.insert(2, "test_observed_pairs", int(test_y.sum()))
    metrics.insert(3, "test_prevalence", float(test_y.mean()))

    metrics.to_csv(args.output_dir / "support_zero_exclusion_metrics.tsv", sep="\t", index=False)
    effects.to_csv(args.output_dir / "support_zero_exclusion_incremental_effects.tsv", sep="\t", index=False)
    bootstrap.to_csv(args.output_dir / "support_zero_exclusion_bootstrap_ap.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(selection_rows).to_csv(args.output_dir / "support_zero_exclusion_hyperparameters.tsv", sep="\t", index=False)
    split_counts.to_csv(args.output_dir / "support_zero_exclusion_counts.tsv", sep="\t", index=False)

    report = {
        "status": "pass",
        "analysis_status": "post hoc sensitivity analysis developed in response to external feedback",
        "exclusion": "other_context_support == 0 removed from train, validation, and test",
        "full_pairs": int(len(full_grid)),
        "retained_pairs": int(len(reduced)),
        "excluded_pairs": int(len(full_grid) - len(reduced)),
        "test_pairs": int(test.sum()),
        "test_observed_pairs": int(test_y.sum()),
        "test_prevalence": float(test_y.mean()),
        "test_components": int(predictions["cluster_representative"].nunique()),
        "bootstrap_replicates": args.bootstrap,
        "seed": args.seed,
    }
    (args.output_dir / "support_zero_exclusion_audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print(effects.to_string(index=False))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
