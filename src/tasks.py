import asyncio

from fastapi.encoders import jsonable_encoder

from configs.database import get_db_session
from configs.redis import set_cache
from src.celery_worker import celery_app
from src.dao.molecule_dao import MoleculeDAO


@celery_app.task
def fetch_and_cache_molecules(identifier: str, limit: int):
    async def run_async():
        async for session in get_db_session():
            dao = MoleculeDAO(session=session)
            cache_key = f"molecules:{identifier}:{limit}"

            try:
                molecules = [
                    molecule async for molecule in
                    dao.list_molecules(identifier=identifier, limit=limit)
                ]
                set_cache(cache_key, molecules)
                return jsonable_encoder(molecules),
            except Exception as e:
                return {"status": "error", "error": str(e)}

    return asyncio.run(run_async())
