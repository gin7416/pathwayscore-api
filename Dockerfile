FROM python:3.11-slim
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .
EXPOSE 8000
# Local/research demo; do not expose publicly with identifiable patient data.
CMD ["uvicorn", "pathwayscore_service.app:app", "--host", "0.0.0.0", "--port", "8000"]
