"""FastAPI research-model serving demo; local/offline by default."""
from typing import Any
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, ConfigDict, Field
from .model import InputValidationError, catalog, canonical_cancer, requirements, score

WARNING = "Research/engineering demonstration only. Not for patient care, medical decisions, or absolute risk prediction."

app = FastAPI(
    title="PathwayScore API",
    description="Transparent fixed-formula Cox scoring on precomputed pathway features; " + WARNING,
    version="0.1.0",
)


class ScoreRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    cancer: str = Field(description="IDC (alias BRCA), LUAD, PAAD, PRAD, or CRC")
    features: dict[str, Any] = Field(description="Use GET /models/{cancer} for feature names and semantics")


class BatchRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    samples: list[ScoreRequest] = Field(min_length=1, max_length=100)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "pathwayscore-api", "version": "0.1.0"}


@app.get("/models")
def models() -> dict[str, Any]:
    return {"cancers": list(catalog()["cancers"]), "alias": {"IDC": "BRCA"}, "warning": WARNING}


@app.get("/models/{cancer}")
def model_details(cancer: str) -> dict[str, Any]:
    try:
        return requirements(cancer)
    except InputValidationError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.post("/score")
def single_score(request: ScoreRequest) -> dict[str, Any]:
    try:
        return score(request.cancer, request.features)
    except InputValidationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.post("/batch")
def batch_score(request: BatchRequest) -> dict[str, Any]:
    results = []
    for index, item in enumerate(request.samples):
        try:
            results.append(score(item.cancer, item.features))
        except InputValidationError as exc:
            raise HTTPException(status_code=422, detail=f"samples[{index}]: {exc}") from exc
    return {"count": len(results), "results": results, "warning": WARNING}
