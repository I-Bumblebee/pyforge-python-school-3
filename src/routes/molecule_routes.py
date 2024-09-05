import logging
from typing import Optional

from fastapi import APIRouter, Depends, Query
from fastapi.encoders import jsonable_encoder
from pydantic import BaseModel

from configs.redis import get_cached_result, set_cache
from configs.setup_logger import setup_logger
from dao.molecule_dao import MoleculeDAO

logger = logging.getLogger(__name__)
setup_logger()

router = APIRouter()


class MoleculeCreate(BaseModel):
    identifier: str
    smiles: str


@router.post("/molecules/")
async def store_molecule(
    molecule: MoleculeCreate, dao: MoleculeDAO = Depends()
):
    logger.info(
        "Received request to create molecule",
        extra={"identifier": molecule.identifier}
    )
    res = await dao.create_molecule(molecule.identifier, molecule.smiles)
    logger.info(
        "Molecule created successfully",
        extra={"identifier": molecule.identifier}
    )
    return res


@router.get("/molecules/{identifier}")
async def view_molecule(identifier: str, dao: MoleculeDAO = Depends()):
    logger.info(
        "Received request to view molecule", extra={"identifier": identifier}
    )

    cache_key = f"molecule:{identifier}"

    cached_result = get_cached_result(cache_key)
    if cached_result:
        logger.info("Cache hit for molecule", extra={"identifier": identifier})
        return cached_result

    try:
        res = await dao.get_molecule_by_identifier(identifier)

        logger.info(
            "Molecule retrieved successfully",
            extra={"identifier": identifier}
        )

        set_cache(cache_key, jsonable_encoder(res))

    except Exception as e:
        logger.error(
            "Error retrieving molecule",
            exc_info=e,
            extra={"identifier": identifier}
        )
        raise

    return res


@router.put("/molecules/{identifier}")
async def update_molecule(
    identifier: str, molecule: MoleculeCreate, dao: MoleculeDAO = Depends()
):
    logger.info(
        "Received request to update molecule",
        extra={"identifier": identifier}
    )
    try:
        res = await dao.update_molecule(identifier, molecule.smiles)
        logger.info(
            "Molecule updated successfully", extra={"identifier": identifier}
        )
    except Exception as e:
        logger.error(
            "Error updating molecule",
            exc_info=e,
            extra={"identifier": identifier}
        )
        raise
    return res


@router.delete("/molecules/{identifier}")
async def destroy_molecule(identifier: str, dao: MoleculeDAO = Depends()):
    logger.info(
        "Received request to delete molecule",
        extra={"identifier": identifier}
    )
    try:
        await dao.delete_molecule(identifier)
        logger.info(
            "Molecule deleted successfully", extra={"identifier": identifier}
        )
    except Exception as e:
        logger.error(
            "Error deleting molecule",
            exc_info=e,
            extra={"identifier": identifier}
        )
        raise
    return {"detail": "Molecule deleted"}


@router.get("/molecules/")
async def index_molecules(
    identifier: Optional[str] = None,
    limit: Optional[int] = Query(None, le=1000),
    dao: MoleculeDAO = Depends()
):
    logger.info(
        "Received request to list molecules",
        extra={
            "identifier": identifier,
            "limit": limit
        }
    )

    cache_key = f"molecules:{identifier}:{limit}"

    cached_result = get_cached_result(cache_key)
    if cached_result is not None:
        logger.info("Cache hit for molecules", extra={"cache_key": cache_key})
        return cached_result

    try:
        molecules = [
            molecule async for molecule in
            dao.list_molecules(identifier=identifier, limit=limit)
        ]
        logger.info(
            "Molecules listed successfully", extra={"count": len(molecules)}
        )

        set_cache(cache_key, jsonable_encoder(molecules))

    except Exception as e:
        logger.error(
            "Error listing molecules",
            exc_info=e,
            extra={
                "identifier": identifier,
                "limit": limit
            }
        )
        raise

    return molecules


@router.get("/molecules/search/")
async def search_molecules_by_substructure(
    identifier: Optional[str] = Query(
        None,
        description=(
            "Identifier of a molecule whose substructure will be used to "
            "find and match other molecules with similar substructures."
        )
    ),
    dao: MoleculeDAO = Depends()
):
    logger.info(
        "Received request to search molecules by substructure",
        extra={"identifier": identifier}
    )
    try:
        res = await dao.find_molecules_by_substructure(identifier)
        logger.info(
            "Molecules found by substructure search",
            extra={
                "identifier": identifier,
                "count": len(res)
            }
        )
    except Exception as e:
        logger.error(
            "Error searching molecules by substructure",
            exc_info=e,
            extra={"identifier": identifier}
        )
        raise
    return res
