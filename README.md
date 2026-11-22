# PathwayScore API — from research model to testable software service

[![tests](https://github.com/gin7416/pathwayscore-api/actions/workflows/ci.yml/badge.svg)](https://github.com/gin7416/pathwayscore-api/actions/workflows/ci.yml)

[繁體中文](#中文說明) | [English](#english-overview)

## 中文說明

**專案定位：** 這是以碩士研究 *Multi-variant pathway aggregation of sparse targeted-panel alterations for prognostic modeling in five solid tumors* 為基礎的**軟體工程延伸展示**，將既有、固定係數的 Cox pathway risk score 包裝成 REST API；**不是重新訓練模型，也不是臨床使用工具**。

原始研究程式碼與固定公式：
https://github.com/ICLab-pathway-level-risk-scores/pathway-level-risk-scores

### 專案動機

研究程式原本側重批次分析與論文結果重現。這個延伸讓研究模型具備：

1. 可查詢的模型輸入規格（`GET /models/{cancer}`）
2. 具輸入驗證的單筆與批次模型評分（`POST /score`, `POST /batch`）
3. 固定的研究模型係數與訓練集標準化參數（`src/.../models.json`），並以自動化測試比對上游研究設定
4. 自動化測試（pytest、GitHub Actions，Python 3.11–3.13）
5. 可重現執行環境（Docker、Docker Compose，CI 中實際建置與啟動）

**僅接受已完成 pathway 特徵工程的數值**。基因突變、CNA、SV 到 pathway 特徵的轉換仍由上游原始研究程式 `preprocessing/build_pathway_features.py` 完成。本服務沒有聲稱可以直接讀原始基因變異檔。

### 快速啟動（不需要 Docker）

在 Python 3.11+ 環境，於專案資料夾：

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
uvicorn pathwayscore_service.app:app --reload
```

瀏覽 http://127.0.0.1:8000/docs 開啟互動 API 文件。

### Docker 啟動

```bash
docker compose up --build
```

API 文件：http://127.0.0.1:8000/docs

> 本展示專案只將 Docker 綁定至本機 `127.0.0.1`，不得使用真實病患識別資訊，**未設定正式對外服務的驗證、權限與安全機制**。

### API 使用範例

先看輸入欄位：

```bash
curl http://127.0.0.1:8000/models/LUAD
```

再用**完全虛構、非真實病患資料**做示範：

```bash
curl -X POST http://127.0.0.1:8000/score \
  -H "Content-Type: application/json" \
  --data-binary @examples/luad_synthetic.json
```

輸出 `log_score`、`relative_score` 及每項特徵貢獻，**沒有預測存活月份或個人的絕對存活率**。

### 特徵工程與科學限制

- 模型使用原研究的五種癌症固定係數：BRCA（論文稱 IDC）、LUAD、PAAD、PRAD、CRC。
- 對 `*_rate_z` 類模型特徵，API 的輸入是相對應 `*_rate_raw`（0~1），由 API 使用原研究 MSK 訓練集的 `mean/std` 做固定標準化。
- `*_zsum` 必須是**已由上游研究定義算好的 composite 值**，不等於單純變異數；`PW_HIT_COUNT` 也必須依上游 pathway 定義計算。API 不替使用者推斷基因層級資料。
- `relative_score` 是 `exp(sum(beta * x))` 的 Cox partial hazard：可在一致處理的研究資料中比較模型分數，但**不是**病患的生存機率、診斷結果或治療建議。
- 原論文顯示外部驗證時不同癌症的 discrimination 與 calibration 表現不同，特別是 CRC 外部結果較弱；本專案不提供跨機構通用的臨床高低風險標籤。
- 本專案的係數及正規化數值來自上游公開 `compact_score_formula.yaml`。`tests/test_upstream_parity.py` 會將 `models.json` 與固定版本（commit `2ef12c6`）的上游 YAML 自動比對：5 種癌症、24 項特徵的係數、轉換類型與標準化參數須**完全一致**。測試驗證的是此服務的數值計算、與研究設定的一致性、API 行為及錯誤處理，**不是**獨立的臨床效能重現。

### 自動化驗證

```bash
pytest -q --cov=pathwayscore_service
```

測試包含五種癌症的特徵介面、獨立公式數值比對、與上游研究設定的參數一致性、別名、缺漏欄位、未知欄位、非數字、無限大／NaN、邊界、批次回應、Swagger 文件等。

GitHub Actions（`.github/workflows/ci.yml`）在每次 push／PR 執行兩個工作：

- **tests**：於 Python 3.11、3.12、3.13 執行完整 pytest 並輸出覆蓋率。
- **docker**：以 Docker Compose 建置並啟動容器，確認 `/health` 正常、`/score` 對範例資料回傳預期分數、缺漏欄位回傳 HTTP 422。

`results/local_pytest_report.txt` 為本機（macOS）執行紀錄。

## English overview

A reproducible engineering extension of a fixed-coefficient pathway-level Cox research model. Provides transparent scoring, feature validation, batch inference, FastAPI docs, pytest coverage and Docker packaging. Inputs are **precomputed pathway features**, not raw genomics; the API is **not for clinical decision-making**. Research coefficients and MSK-derived standardization parameters come from the MIT-licensed original model configuration and are checked against a pinned copy of the upstream YAML by an automated parity test. CI runs the test suite on Python 3.11–3.13 and builds, starts and smoke-tests the Docker container. Refer to the upstream research repository for the original scientific analysis and data access.

### Source and attribution

Original manuscript and model repository: https://github.com/ICLab-pathway-level-risk-scores/pathway-level-risk-scores

Original code license: MIT, Copyright (c) 2026 Kai-Ling Kao and Shinn-Ying Ho. The original code and this separate engineering extension should not be represented as identical repositories. No patient-level research data is included in this repository. `tests/fixtures/upstream_compact_score_formula.yaml` is an unmodified copy of the upstream formula file (commit `2ef12c6d445c1598ee5ffb4f0885ba91e87308ca`), included under its MIT license for parity testing.
