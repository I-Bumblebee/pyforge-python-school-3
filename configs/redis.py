import redis
import json
from .settings import settings

redis_client = redis.Redis(port=settings.redis_port)


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(
    key: str, value: dict | list, expiration: int = settings.cache_duration
):
    redis_client.setex(key, expiration, json.dumps(value))
