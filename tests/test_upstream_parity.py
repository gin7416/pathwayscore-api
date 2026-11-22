"""Parity check: service model config must match the upstream research formula.

`tests/fixtures/upstream_compact_score_formula.yaml` is an unmodified copy of
compact_score/compact_score_formula.yaml from
https://github.com/ICLab-pathway-level-risk-scores/pathway-level-risk-scores
pinned at commit 2ef12c6d445c1598ee5ffb4f0885ba91e87308ca (MIT License).
"""
from pathlib import Path

import pytest
import yaml

from pathwayscore_service.model import catalog

UPSTREAM_PATH = Path(__file__).parent / "fixtures" / "upstream_compact_score_formula.yaml"
UPSTREAM = yaml.safe_load(UPSTREAM_PATH.read_text(encoding="utf-8"))
SERVICE = catalog()


def test_same_cancers_and_alias():
    assert set(SERVICE["cancers"]) == set(UPSTREAM["cancers"])
    assert SERVICE["aliases"] == UPSTREAM["cancer_key_aliases"]


@pytest.mark.parametrize("cancer", sorted(UPSTREAM["cancers"]))
def test_feature_set_and_pathways_match(cancer):
    upstream, service = UPSTREAM["cancers"][cancer], SERVICE["cancers"][cancer]
    assert service["pathways"] == upstream["hit_count_pathways"]
    assert [f["name"] for f in service["features"]] == [f["name"] for f in upstream["features"]]
    assert len(service["features"]) == upstream["K"]


@pytest.mark.parametrize("cancer", sorted(UPSTREAM["cancers"]))
def test_coefficients_and_standardization_match_exactly(cancer):
    pairs = zip(UPSTREAM["cancers"][cancer]["features"], SERVICE["cancers"][cancer]["features"])
    for upstream, service in pairs:
        assert service["beta"] == upstream["beta"], upstream["name"]
        assert service["transform"] == upstream["transform"], upstream["name"]
        assert service.get("mean") == upstream.get("zscore_mean"), upstream["name"]
        assert service.get("std") == upstream.get("zscore_std"), upstream["name"]
