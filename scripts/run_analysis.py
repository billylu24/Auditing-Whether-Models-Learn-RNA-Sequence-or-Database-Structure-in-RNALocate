"""Run the scoped RNALocate context-shortcut analyses for the high-school paper.

The outcome in the primary analysis is database record availability, not
biological localization. All new outputs are isolated under results/high_school.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]

MODEL_FEATURES: dict[str, tuple[str, ...]] = {
    "prevalence": (),
    "context_only": ("context",),
    "support_only": ("support",),
    "sequence_only": ("sequence",),
    "context_sequence": ("context", "sequence"),
    "context_support": ("context", "support"),
    "support_sequence": ("support", "sequence"),
    "full": ("context", "support", "sequence"),
}

DISPLAY_NAMES = {
    "prevalence": "Prevalence",
    "context_only": "Context",
    "support_only": "Other-context support",
    "sequence_only": "Low-level sequence",
    "context_sequence": "Context + sequence",
    "context_support": "Context + support",
    "support_sequence": "Support + sequence",
    "full": "Context + support + sequence",
}

FEATURE_VIEWS: dict[str, np.ndarray] = {
    "composition_length": np.arange(64, 70),
    "three_mer_only": np.arange(0, 64),
    "complete_70": np.arange(0, 70),
}

CS = (0.01, 0.1, 1.0, 10.0)
CONTRASTS = {
    "sequence_after_context": ("context_sequence", "context_only"),
    "support_after_context": ("context_support", "context_only"),
    "sequence_after_context_support": ("full", "context_support"),
    "support_after_context_sequence": ("full", "context_sequence"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        type=Path,
        default=ROOT / "inputs/observed_gene_context_pairs.tsv.gz",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=ROOT / "inputs/eligible_genes.tsv.gz",
    )
    parser.add_argument(
        "--split",
        type=Path,
        default=ROOT / "inputs/gene_component_split.tsv",
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=ROOT / "inputs/sequence_features_70.npy",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--calibration-bins", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260814)
    return parser.parse_args()


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return math.nan, math.nan
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    radius = z * math.sqrt(proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def safe_metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    result = {
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }
    result["roc_auc"] = float(roc_auc_score(y, probability)) if np.unique(y).size == 2 else math.nan
    return result


def build_grid(args: argparse.Namespace) -> tuple[pd.DataFrame, np.ndarray]:
    observed = pd.read_csv(args.panel, sep="\t", usecols=["rna_symbol", "context"]).drop_duplicates()
    reference = pd.read_csv(args.reference, sep="\t", usecols=["rna_symbol"]).drop_duplicates()
    genes = reference["rna_symbol"].tolist()
    contexts = sorted(observed["context"].unique())
    cached = np.load(args.features, mmap_mode="r")
    if cached.shape != (len(genes), 70):
        raise ValueError(f"Expected a ({len(genes)}, 70) sequence cache, found {cached.shape}")
    if set(observed["rna_symbol"]) != set(genes):
        raise ValueError("Observation panel and feature-reference gene sets differ")

    grid = pd.MultiIndex.from_product([genes, contexts], names=["rna_symbol", "context"]).to_frame(index=False)
    grid = grid.merge(observed.assign(observed=1), on=["rna_symbol", "context"], how="left", validate="one_to_one")
    grid["observed"] = grid["observed"].fillna(0).astype(np.uint8)
    total_support = observed.groupby("rna_symbol", observed=True)["context"].nunique()
    grid["other_context_support"] = (
        grid["rna_symbol"].map(total_support).astype(np.int16) - grid["observed"].astype(np.int16)
    )

    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    context_lookup = {context: index for index, context in enumerate(contexts)}
    grid["gene_index"] = grid["rna_symbol"].map(gene_lookup).astype(np.int32)
    grid["context_index"] = grid["context"].map(context_lookup).astype(np.int16)
    split = pd.read_csv(args.split, sep="\t", usecols=["rna_symbol", "cluster_representative", "split"])
    grid = grid.merge(split, on="rna_symbol", how="left", validate="many_to_one")
    if grid[["cluster_representative", "split"]].isna().any().any():
        raise ValueError("Missing split assignment on the complete grid")
    return grid, cached


def design_matrix(
    parts: tuple[str, ...],
    gene_index: np.ndarray,
    context_index: np.ndarray,
    support: np.ndarray,
    sequence_features: np.ndarray,
    context_count: int,
    sequence_columns: np.ndarray | None = None,
) -> np.ndarray:
    values: list[np.ndarray] = []
    if "context" in parts:
        values.append(np.eye(context_count, dtype=np.float32)[context_index])
    if "support" in parts:
        values.append(support[:, None].astype(np.float32))
    if "sequence" in parts:
        columns = FEATURE_VIEWS["complete_70"] if sequence_columns is None else sequence_columns
        values.append(np.asarray(sequence_features[gene_index][:, columns], dtype=np.float32))
    if not values:
        return np.empty((len(gene_index), 0), dtype=np.float32)
    return np.concatenate(values, axis=1)


@dataclass
class FittedResult:
    probability: np.ndarray
    selected_c: float | None
    validation_ap: float | None
    selection_rows: list[dict[str, object]]


def fit_model(
    x: np.ndarray,
    y: np.ndarray,
    train: np.ndarray,
    validation: np.ndarray,
    test: np.ndarray,
    model_name: str,
    seed: int,
) -> FittedResult:
    fit = train | validation
    if x.shape[1] == 0:
        probability = np.full(int(test.sum()), float(y[fit].mean()), dtype=np.float64)
        return FittedResult(probability, None, None, [])

    selection_rows: list[dict[str, object]] = []
    best_c: float | None = None
    best_ap = -np.inf
    for c_value in CS:
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=c_value, solver="lbfgs", max_iter=750, random_state=seed),
        )
        model.fit(x[train], y[train])
        probability = model.predict_proba(x[validation])[:, 1]
        ap = float(average_precision_score(y[validation], probability))
        selection_rows.append(
            {
                "model": model_name,
                "C": c_value,
                "validation_average_precision": ap,
                "selected": False,
            }
        )
        if ap > best_ap:
            best_c, best_ap = c_value, ap
    assert best_c is not None
    for row in selection_rows:
        if row["C"] == best_c:
            row["selected"] = True

    final_model = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=best_c, solver="lbfgs", max_iter=750, random_state=seed),
    )
    final_model.fit(x[fit], y[fit])
    probability = final_model.predict_proba(x[test])[:, 1]
    return FittedResult(probability, best_c, best_ap, selection_rows)


def fit_factorial_models(
    grid: pd.DataFrame,
    cached: np.ndarray,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    train = grid["split"].eq("train").to_numpy()
    validation = grid["split"].eq("validation").to_numpy()
    test = grid["split"].eq("test").to_numpy()
    y = grid["observed"].to_numpy(dtype=np.uint8)
    gene = grid["gene_index"].to_numpy()
    context = grid["context_index"].to_numpy()
    support = grid["other_context_support"].to_numpy()
    context_count = int(context.max()) + 1

    predictions = grid.loc[
        test,
        ["rna_symbol", "context", "cluster_representative", "observed", "other_context_support"],
    ].reset_index(drop=True)
    selection_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []

    for model_name, parts in MODEL_FEATURES.items():
        started = time.perf_counter()
        x = design_matrix(parts, gene, context, support, cached, context_count)
        fitted = fit_model(x, y, train, validation, test, model_name, seed)
        predictions[f"prob_{model_name}"] = fitted.probability
        selection_rows.extend(fitted.selection_rows)
        values = safe_metrics(y[test], fitted.probability)
        metric_rows.append(
            {
                "model": model_name,
                "display_name": DISPLAY_NAMES[model_name],
                "selected_C": fitted.selected_c,
                "validation_average_precision": fitted.validation_ap,
                **values,
                "fit_seconds": time.perf_counter() - started,
            }
        )
        del x
        gc.collect()
        print(f"factorial {model_name}: AP={values['average_precision']:.6f}", flush=True)

    return predictions, pd.DataFrame(metric_rows), pd.DataFrame(selection_rows)


class WeightedAveragePrecision:
    def __init__(self, y: np.ndarray, probability: np.ndarray):
        self.order = np.argsort(-probability, kind="mergesort")
        sorted_probability = probability[self.order]
        self.sorted_y = y[self.order].astype(np.float64)
        self.starts = np.r_[0, np.flatnonzero(sorted_probability[1:] != sorted_probability[:-1]) + 1]

    def score(self, row_weights: np.ndarray) -> float:
        weights = row_weights[self.order].astype(np.float64, copy=False)
        group_total = np.add.reduceat(weights, self.starts)
        group_positive = np.add.reduceat(weights * self.sorted_y, self.starts)
        total_positive = group_positive.sum()
        if total_positive <= 0:
            return math.nan
        cumulative_total = np.cumsum(group_total)
        cumulative_positive = np.cumsum(group_positive)
        precision = np.divide(
            cumulative_positive,
            cumulative_total,
            out=np.zeros_like(cumulative_positive),
            where=cumulative_total > 0,
        )
        return float(np.sum(precision * group_positive) / total_positive)


def bootstrap_average_precision(
    predictions: pd.DataFrame,
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    components, component_index = np.unique(predictions["cluster_representative"].to_numpy(), return_inverse=True)
    y = predictions["observed"].to_numpy(dtype=np.uint8)
    calculators = {
        model: WeightedAveragePrecision(y, predictions[f"prob_{model}"].to_numpy())
        for model in MODEL_FEATURES
    }
    values = np.empty((replicates, len(MODEL_FEATURES)), dtype=np.float64)
    rng = np.random.default_rng(seed)
    for replicate in range(replicates):
        sampled = rng.integers(0, len(components), size=len(components))
        component_weights = np.bincount(sampled, minlength=len(components))
        row_weights = component_weights[component_index]
        for column, calculator in enumerate(calculators.values()):
            values[replicate, column] = calculator.score(row_weights)
        if (replicate + 1) % 100 == 0:
            print(f"bootstrap {replicate + 1}/{replicates}", flush=True)

    model_names = list(MODEL_FEATURES)
    bootstrap = pd.DataFrame(values, columns=model_names)
    interval_rows = []
    for model in model_names:
        interval_rows.append(
            {
                "model": model,
                "ap_ci95_low": float(np.nanquantile(bootstrap[model], 0.025)),
                "ap_ci95_high": float(np.nanquantile(bootstrap[model], 0.975)),
                "bootstrap_replicates": replicates,
                "bootstrap_unit": "all-sequence graph component",
            }
        )

    effect_rows = []
    for contrast, (left, right) in CONTRASTS.items():
        difference = bootstrap[left] - bootstrap[right]
        effect_rows.append(
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
    return bootstrap, pd.DataFrame(interval_rows), pd.DataFrame(effect_rows)


def coverage_gradient(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pooled_rows: list[dict[str, object]] = []
    context_rows: list[dict[str, object]] = []
    for support, frame in predictions.groupby("other_context_support", sort=True):
        successes = int(frame["observed"].sum())
        total = int(len(frame))
        low, high = wilson_interval(successes, total)
        pooled_rows.append(
            {
                "other_context_support": int(support),
                "test_pairs": total,
                "observed_pairs": successes,
                "observed_fraction": successes / total,
                "ci95_low": low,
                "ci95_high": high,
                "selection_boundary": int(support) == 0,
            }
        )
    for (context, support), frame in predictions.groupby(["context", "other_context_support"], sort=True):
        successes = int(frame["observed"].sum())
        total = int(len(frame))
        low, high = wilson_interval(successes, total)
        context_rows.append(
            {
                "context": context,
                "other_context_support": int(support),
                "test_pairs": total,
                "observed_pairs": successes,
                "observed_fraction": successes / total,
                "ci95_low": low,
                "ci95_high": high,
            }
        )
    return pd.DataFrame(pooled_rows), pd.DataFrame(context_rows)


def permute_support_within_context_and_split(
    support: np.ndarray,
    context: np.ndarray,
    split_codes: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    shuffled = support.copy()
    for split_code in np.unique(split_codes):
        for context_code in np.unique(context):
            local = np.flatnonzero((split_codes == split_code) & (context == context_code))
            shuffled[local] = support[rng.permutation(local)]
    return shuffled


def shuffled_support_control(
    grid: pd.DataFrame,
    real_ap: float,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    train = grid["split"].eq("train").to_numpy()
    validation = grid["split"].eq("validation").to_numpy()
    test = grid["split"].eq("test").to_numpy()
    y = grid["observed"].to_numpy(dtype=np.uint8)
    context = grid["context_index"].to_numpy()
    support = grid["other_context_support"].to_numpy()
    split_lookup = {name: code for code, name in enumerate(sorted(grid["split"].unique()))}
    split_codes = grid["split"].map(split_lookup).to_numpy(dtype=np.int8)
    context_one_hot = np.eye(int(context.max()) + 1, dtype=np.float32)[context]
    rows = []

    for replicate in range(replicates):
        rng = np.random.default_rng(seed + replicate)
        shuffled = permute_support_within_context_and_split(support, context, split_codes, rng)
        x = np.column_stack([context_one_hot, shuffled.astype(np.float32)])
        fitted = fit_model(x, y, train, validation, test, f"shuffled_support_{replicate:03d}", seed + replicate)
        result = safe_metrics(y[test], fitted.probability)
        rows.append(
            {
                "replicate": replicate,
                "seed": seed + replicate,
                "selected_C": fitted.selected_c,
                "validation_average_precision": fitted.validation_ap,
                **result,
                "real_context_support_ap": real_ap,
                "real_minus_shuffled_ap": real_ap - result["average_precision"],
            }
        )
        if (replicate + 1) % 10 == 0:
            print(f"support permutation {replicate + 1}/{replicates}", flush=True)
    return pd.DataFrame(rows)


def context_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    models = ("prevalence", "context_only", "context_support", "full")
    for context, frame in predictions.groupby("context", sort=True):
        y = frame["observed"].to_numpy(dtype=np.uint8)
        prevalence = float(y.mean())
        for model in models:
            probability = frame[f"prob_{model}"].to_numpy()
            result = safe_metrics(y, probability)
            denominator = 1.0 - prevalence
            rows.append(
                {
                    "context": context,
                    "model": model,
                    "display_name": DISPLAY_NAMES[model],
                    "test_pairs": int(len(frame)),
                    "observed_pairs": int(y.sum()),
                    "prevalence": prevalence,
                    **result,
                    "ap_lift": result["average_precision"] - prevalence,
                    "normalized_ap_lift": (
                        (result["average_precision"] - prevalence) / denominator if denominator > 0 else math.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def calibration_bins(predictions: pd.DataFrame, bins: int) -> pd.DataFrame:
    rows = []
    y = predictions["observed"].to_numpy(dtype=np.uint8)
    for model in ("context_only", "context_support", "full"):
        probability = predictions[f"prob_{model}"].to_numpy()
        edges = np.unique(np.quantile(probability, np.linspace(0.0, 1.0, bins + 1)))
        if len(edges) == 1:
            bin_index = np.zeros(len(probability), dtype=np.int16)
        else:
            bin_index = np.searchsorted(edges[1:-1], probability, side="right")
        for local_bin in np.unique(bin_index):
            take = bin_index == local_bin
            successes = int(y[take].sum())
            total = int(take.sum())
            low, high = wilson_interval(successes, total)
            rows.append(
                {
                    "model": model,
                    "display_name": DISPLAY_NAMES[model],
                    "bin": int(local_bin) + 1,
                    "test_pairs": total,
                    "mean_predicted_probability": float(probability[take].mean()),
                    "observed_fraction": successes / total,
                    "ci95_low": low,
                    "ci95_high": high,
                }
            )
    return pd.DataFrame(rows)


def sequence_sensitivity(
    grid: pd.DataFrame,
    cached: np.ndarray,
    seed: int,
) -> pd.DataFrame:
    train = grid["split"].eq("train").to_numpy()
    validation = grid["split"].eq("validation").to_numpy()
    test = grid["split"].eq("test").to_numpy()
    y = grid["observed"].to_numpy(dtype=np.uint8)
    gene = grid["gene_index"].to_numpy()
    context = grid["context_index"].to_numpy()
    support = grid["other_context_support"].to_numpy()
    context_count = int(context.max()) + 1
    rows = []
    for view_name, columns in FEATURE_VIEWS.items():
        for model_name, parts in (("sequence_only", ("sequence",)), ("full", ("context", "support", "sequence"))):
            started = time.perf_counter()
            x = design_matrix(parts, gene, context, support, cached, context_count, sequence_columns=columns)
            fitted = fit_model(x, y, train, validation, test, f"{model_name}_{view_name}", seed)
            result = safe_metrics(y[test], fitted.probability)
            rows.append(
                {
                    "feature_view": view_name,
                    "feature_count": int(len(columns)),
                    "model": model_name,
                    "selected_C": fitted.selected_c,
                    "validation_average_precision": fitted.validation_ap,
                    **result,
                    "fit_seconds": time.perf_counter() - started,
                }
            )
            del x
            gc.collect()
            print(f"sensitivity {model_name}/{view_name}: AP={result['average_precision']:.6f}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    grid, cached = build_grid(args)

    predictions, metrics, selection = fit_factorial_models(grid, cached, args.seed)
    bootstrap, intervals, effects = bootstrap_average_precision(predictions, args.bootstrap, args.seed + 1000)
    metrics = metrics.merge(intervals, on="model", how="left", validate="one_to_one")
    metric_lookup = metrics.set_index("model")["average_precision"]
    effect_records = effects.to_dict("records")
    for row in effect_records:
        left = row["left_model"]
        right = row["right_model"]
        row["ap_difference"] = float(metric_lookup.loc[left] - metric_lookup.loc[right])
    effects = pd.DataFrame(effect_records)

    pooled_gradient, context_gradient = coverage_gradient(predictions)
    permutation = shuffled_support_control(
        grid,
        float(metrics.set_index("model").loc["context_support", "average_precision"]),
        args.permutations,
        args.seed + 2000,
    )
    per_context = context_metrics(predictions)
    calibration = calibration_bins(predictions, args.calibration_bins)
    sensitivity = sequence_sensitivity(grid, cached, args.seed)

    predictions.to_csv(args.output_dir / "record_availability_test_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    metrics.to_csv(args.output_dir / "record_availability_metrics.tsv", sep="\t", index=False)
    selection.to_csv(args.output_dir / "record_availability_hyperparameters.tsv", sep="\t", index=False)
    bootstrap.to_csv(args.output_dir / "record_availability_bootstrap_ap.tsv.gz", sep="\t", index=False, compression="gzip")
    effects.to_csv(args.output_dir / "record_availability_incremental_effects.tsv", sep="\t", index=False)
    pooled_gradient.to_csv(args.output_dir / "support_gradient.tsv", sep="\t", index=False)
    context_gradient.to_csv(args.output_dir / "support_gradient_by_context.tsv", sep="\t", index=False)
    permutation.to_csv(args.output_dir / "shuffled_support_control.tsv", sep="\t", index=False)
    per_context.to_csv(args.output_dir / "record_availability_by_context.tsv", sep="\t", index=False)
    calibration.to_csv(args.output_dir / "calibration_bins.tsv", sep="\t", index=False)
    sensitivity.to_csv(args.output_dir / "sequence_view_sensitivity.tsv", sep="\t", index=False)

    support_zero = pooled_gradient.loc[pooled_gradient["other_context_support"].eq(0)]
    permutation_p = float(
        (1 + (permutation["average_precision"] >= permutation["real_context_support_ap"]).sum())
        / (1 + len(permutation))
    )
    report = {
        "status": "pass",
        "analysis_scope": "RNALocate record availability and context-shortcut audit",
        "biological_negative_interpretation": False,
        "grid_rows": int(len(grid)),
        "genes": int(grid["rna_symbol"].nunique()),
        "contexts": int(grid["context"].nunique()),
        "observed_pairs": int(grid["observed"].sum()),
        "unobserved_pairs": int((1 - grid["observed"]).sum()),
        "test_rows": int(len(predictions)),
        "test_components": int(predictions["cluster_representative"].nunique()),
        "factorial_models": list(MODEL_FEATURES),
        "bootstrap_replicates": args.bootstrap,
        "permutation_replicates": args.permutations,
        "permutation_empirical_p_upper_tail": permutation_p,
        "support_zero_is_all_observed": bool(
            len(support_zero) == 1 and float(support_zero.iloc[0]["observed_fraction"]) == 1.0
        ),
        "support_zero_explanation": (
            "Eligible genes have at least one retained record; zero support in all other contexts "
            "therefore implies a record in the candidate context."
        ),
        "runtime_seconds": time.perf_counter() - started,
        "seed": args.seed,
    }
    (args.output_dir / "audit_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
