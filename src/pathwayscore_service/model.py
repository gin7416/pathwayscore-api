"""Fixed-coefficient pathway-level scoring with strict input validation.

This is a software engineering extension to the cited research repository, not a
replacement for upstream pathway feature engineering and not a clinical tool.
"""
import json
import math
from functools import lru_cache
from importlib.resources import files
from typing import Any


class InputValidationError(ValueError):
    """Invalid model inputs; safe to return as HTTP 422."""


@lru_cache(maxsize=1)
def catalog() -> dict[str, Any]:
    path = files("pathwayscore_service").joinpath("data/models.json")
    result = json.loads(path.read_text(encoding="utf-8"))
    for model in result["cancers"].values():
        assert len(model["features"]) in (4, 5)
    return result


def canonical_cancer(cancer: str) -> str:
    if not isinstance(cancer, str):
        raise InputValidationError("cancer must be a string")
    name = cancer.strip().upper()
    config = catalog()
    name = config["aliases"].get(name, name)
    if name not in config["cancers"]:
        raise InputValidationError(f"Unsupported cancer: {cancer}. Use BRCA/IDC, LUAD, PAAD, PRAD, CRC.")
    return name


def input_field_name(feature: dict[str, Any]) -> str:
    name = feature["name"]
    if feature["transform"] == "zscore_msk_train":
        return name.removesuffix("_z") + "_raw"
    return name


def requirements(cancer: str) -> dict[str, Any]:
    key = canonical_cancer(cancer)
    info = catalog()["cancers"][key]
    fields = []
    for feature in info["features"]:
        kind = feature["transform"]
        fields.append({
            "input_field": input_field_name(feature),
            "model_feature": feature["name"],
            "type": kind,
            "description": {
                "integer_count": f"integer 0..{len(info['pathways'])}, from upstream pathway aggregation",
                "binary": "0 or 1, from upstream pathway aggregation",
                "composite_no_zscore": "precomputed composite from upstream pathway aggregation; NOT a raw count",
                "zscore_msk_train": "raw pathway alteration rate (0..1); normalized with frozen MSK train mean/std"
            }[kind],
        })
    return {"cancer": key, "fields": fields, "input_stage": "precomputed pathway-level features"}


def score(cancer: str, features: dict[str, Any]) -> dict[str, Any]:
    key = canonical_cancer(cancer)
    if not isinstance(features, dict):
        raise InputValidationError("features must be an object")
    model = catalog()["cancers"][key]
    expected = {input_field_name(f) for f in model["features"]}
    missing = sorted(expected - features.keys())
    unknown = sorted(features.keys() - expected)
    if missing or unknown:
        raise InputValidationError(f"Feature keys invalid; missing={missing}, unexpected={unknown}. GET /models/{key} for required inputs.")
    log_score = 0.0
    contributions = []
    for feature in model["features"]:
        field = input_field_name(feature)
        raw = features[field]
        if isinstance(raw, bool) or not isinstance(raw, (int, float)) or not math.isfinite(raw):
            raise InputValidationError(f"{field} must be a finite number")
        kind = feature["transform"]
        if kind == "binary" and raw not in (0, 1):
            raise InputValidationError(f"{field} must be 0 or 1")
        if kind == "integer_count" and (not float(raw).is_integer() or not (0 <= raw <= len(model["pathways"]))):
            raise InputValidationError(f"{field} must be an integer from 0 to {len(model['pathways'])}")
        if kind == "zscore_msk_train" and not (0 <= raw <= 1):
            raise InputValidationError(f"{field} must be a raw alteration fraction between 0 and 1")
        if kind == "composite_no_zscore" and not (-50 <= raw <= 50):
            raise InputValidationError(f"{field} composite is out of supported demo range [-50, 50]")
        transformed = (raw - feature["mean"]) / feature["std"] if kind == "zscore_msk_train" else raw
        effect = feature["beta"] * transformed
        log_score += effect
        contributions.append({
            "feature": feature["name"], "source_field": field,
            "transformed_value": round(transformed, 9), "log_score_contribution": round(effect, 9),
        })
    try:
        partial_hazard = math.exp(log_score)
    except OverflowError:
        raise InputValidationError("Computed score overflow; verify upstream features") from None
    if not math.isfinite(partial_hazard):
        raise InputValidationError("Computed score not finite")
    return {
        "cancer": key,
        "log_score": round(log_score, 9),
        "relative_score": round(partial_hazard, 9),
        "contributions": contributions,
        "interpretation": "Relative Cox partial hazard only; NOT an absolute survival probability, diagnosis, or clinical recommendation.",
    }
