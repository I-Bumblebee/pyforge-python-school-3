FROM python:3.12.5-slim

ARG APP_ENV

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/app/src \
    APP_ENV=$APP_ENV

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

CMD ["alembic -x data=true upgrade head && exec uvicorn src.main:app --host 0.0.0.0 --port 8000 --reload"]