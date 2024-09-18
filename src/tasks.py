import asyncio

from fastapi.encoders import jsonable_encoder

from configs.database import get_db_session
from configs.redis import set_cache
from src.celery_worker import celery_app
from src.dao.molecule_dao import MoleculeDAO


@celery_app.task
def fetch_and_cache_molecules(
    identifier: str,
    limit: int,
):
    dao = MoleculeDAO(session=get_db_session())

    async def run_async():
        cache_key = f"molecules:{identifier}:{limit}"

        try:
            molecules = [
                molecule async for molecule in
                dao.list_molecules(identifier=identifier, limit=limit)
            ]
            set_cache(cache_key, jsonable_encoder(molecules))
            return {
                "status": "success",
                "cache_key": cache_key,
                "molecules": molecules,
            }
        except Exception as e:
            print(f"Error fetching molecules: {e}")
            return {"status": "error", "error": str(e)}

    loop = asyncio.get_event_loop()
    return loop.run_until_complete(run_async())
