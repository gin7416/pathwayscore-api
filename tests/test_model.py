import json
import math
from pathlib import Path
import pytest
from pathwayscore_service.model import InputValidationError, canonical_cancer, catalog, requirements, score

EXAMPLE = json.loads((Path(__file__).parents[1] / "examples/luad_synthetic.json").read_text())


def test_luad_matches_independently_calculated_cox_formula():
    inputs = EXAMPLE["features"]
    expected_log = (
        0.08623152852858001 * 3
        + 0.10802577709803124 * ((0.12 - 0.1260680034873583) / 0.1170473406848421)
        - 0.07991265465222017 * ((0.04 - 0.027615518744551) / 0.0622433492586812)
        + 0.11153169994040811 * 0.1
        - 0.07770157107125727 * ((0.05 - 0.0346095763302899) / 0.0294336706930631)
    )
    result = score(**EXAMPLE)
    assert result["log_score"] == pytest.approx(expected_log, abs=1e-8)
    assert result["relative_score"] == pytest.approx(math.exp(expected_log), abs=1e-8)
    assert len(result["contributions"]) == 5


def test_all_five_cancer_models_have_complete_schema_and_score():
    for cancer in catalog()["cancers"]:
        sample = {}
        for spec in requirements(cancer)["fields"]:
            sample[spec["input_field"]] = 0
        result = score(cancer, sample)
        assert result["relative_score"] > 0
        assert len(result["contributions"]) in (4, 5)


def test_idc_alias():
    assert canonical_cancer("idc") == "BRCA"


def test_missing_feature_rejected():
    with pytest.raises(InputValidationError, match="missing="):
        score("LUAD", {"PW_HIT_COUNT": 1})


def test_unexpected_feature_rejected():
    bad = dict(EXAMPLE["features"], NOT_A_FEATURE=123)
    with pytest.raises(InputValidationError, match="unexpected="):
        score("LUAD", bad)


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), "12"])
def test_non_numeric_and_nonfinite_values_rejected(value):
    bad = dict(EXAMPLE["features"], PW_HIT_COUNT=value)
    with pytest.raises(InputValidationError, match="finite number"):
        score("LUAD", bad)


def test_fraction_out_of_range_rejected():
    bad = dict(EXAMPLE["features"], PW_NOTCH_any_rate_raw=1.2)
    with pytest.raises(InputValidationError, match="between 0 and 1"):
        score("LUAD", bad)


def test_hit_count_bound_uses_cancer_specific_pathways():
    bad = dict(EXAMPLE["features"], PW_HIT_COUNT=10)
    with pytest.raises(InputValidationError, match="0 to 9"):
        score("LUAD", bad)


def test_binary_not_boolean_or_fraction():
    sample = {item["input_field"]: 0 for item in requirements("PAAD")["fields"]}
    sample["PW_RTK_RAS_amp_hit"] = 0.4
    with pytest.raises(InputValidationError, match="must be 0 or 1"):
        score("PAAD", sample)


def test_unsupported_cancer_rejected():
    with pytest.raises(InputValidationError, match="Unsupported"):
        score("OTHER", {})
