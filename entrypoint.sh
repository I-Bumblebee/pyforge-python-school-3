#!/bin/sh

# Run database migrations
alembic -x data=true upgrade head

# Start FastAPI server in the background
uvicorn src.main:app --host 0.0.0.0 --port 8000 --reload &

# Start Celery worker
celery -A src.celery_worker worker --loglevel=info
