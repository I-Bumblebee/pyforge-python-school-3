import json

import redis
from fastapi.encoders import jsonable_encoder

from configs.settings import settings

redis_client = redis.Redis(port=settings.redis_port, host=settings.redis_host)


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(
    key: str,
    value,
    expiration: int = settings.cache_duration,
):
    redis_client.setex(key, expiration, json.dumps(jsonable_encoder(value)))
