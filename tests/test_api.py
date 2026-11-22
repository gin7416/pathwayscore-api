import json
from pathlib import Path
from fastapi.testclient import TestClient
from pathwayscore_service.app import app

client = TestClient(app)
EXAMPLE = json.loads((Path(__file__).parents[1] / "examples/luad_synthetic.json").read_text())


def test_health_and_docs():
    assert client.get("/health").json()["status"] == "ok"
    assert client.get("/docs").status_code == 200
    assert client.get("/openapi.json").status_code == 200


def test_models_schema():
    assert len(client.get("/models").json()["cancers"]) == 5
    response = client.get("/models/IDC")
    assert response.status_code == 200
    assert response.json()["cancer"] == "BRCA"


def test_single_api_score():
    response = client.post("/score", json=EXAMPLE)
    assert response.status_code == 200
    result = response.json()
    assert result["relative_score"] > 0
    assert "NOT an absolute" in result["interpretation"]


def test_batch_api_score():
    response = client.post("/batch", json={"samples": [EXAMPLE, EXAMPLE]})
    assert response.status_code == 200
    assert response.json()["count"] == 2


def test_http_validation_and_batch_index():
    example = dict(EXAMPLE)
    example["features"] = {"PW_HIT_COUNT": 1}
    assert client.post("/score", json=example).status_code == 422
    res = client.post("/batch", json={"samples": [EXAMPLE, example]})
    assert res.status_code == 422
    assert "samples[1]" in res.json()["detail"]


def test_batch_cannot_be_empty():
    assert client.post("/batch", json={"samples": []}).status_code == 422


def test_no_unknown_request_fields():
    assert client.post("/score", json={**EXAMPLE, "medical_decision": True}).status_code == 422
