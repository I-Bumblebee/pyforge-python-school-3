from celery import Celery

from configs.settings import settings

redis_url = f'redis://{settings.redis_host}:{settings.redis_port}/0'

celery_app = Celery(
    'tasks',
    broker=redis_url,
    backend=redis_url,
)

celery_app.conf.update(
    task_track_started=True,
    imports=[
        'src.tasks',
    ],
)
