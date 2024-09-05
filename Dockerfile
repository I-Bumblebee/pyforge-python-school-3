FROM python:3.12.5-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/app/src

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY requirements-dev.txt .

RUN if [ "$APP_ENV" != "production" ]; then \
        pip install --no-cache-dir -r requirements-dev.txt; \
    fi

COPY . .

EXPOSE 8000

ENTRYPOINT ["sh", "-c"]

# When testing do not seed data
CMD if [ "$APP_ENV" = "testing" ]; then \
        alembic upgrade head && exec uvicorn src.main:app --host 0.0.0.0 --port 8000 --reload; \
    else \
        alembic -x data=true upgrade head && exec uvicorn src.main:app --host 0.0.0.0 --port 8000 --reload; \
    fi
