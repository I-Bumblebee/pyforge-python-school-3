from typing import Optional

from fastapi import APIRouter, Depends, Query
from pydantic import BaseModel

from dao.molecule_dao import MoleculeDAO

router = APIRouter()


class MoleculeCreate(BaseModel):
    identifier: str
    smiles: str


@router.post("/molecules/")
async def store_molecule(
        molecule: MoleculeCreate,
        dao: MoleculeDAO = Depends()):
    return await dao.create_molecule(molecule.identifier, molecule.smiles)


@router.get("/molecules/{identifier}")
async def view_molecule(identifier: str, dao: MoleculeDAO = Depends()):
    return await dao.get_molecule_by_identifier(identifier)


@router.put("/molecules/{identifier}")
async def update_molecule(
        identifier: str,
        molecule: MoleculeCreate,
        dao: MoleculeDAO = Depends()):
    return await dao.update_molecule(identifier, molecule.smiles)


@router.delete("/molecules/{identifier}")
async def destroy_molecule(identifier: str, dao: MoleculeDAO = Depends()):
    await dao.delete_molecule(identifier)
    return {"detail": "Molecule deleted"}


@router.get("/molecules/")
async def index_molecules(
        identifier: Optional[str] = None,
        dao: MoleculeDAO = Depends()):
    return await dao.list_molecules(identifier)


@router.get("/molecules/search/")
async def search_molecules_by_substructure(
        identifier: Optional[str] = Query(
            None,
            description="Identifier of a molecule whose\
                  substructure will be used to find and\
                      match other molecules with similar substructures."
        ),
        dao: MoleculeDAO = Depends()):
    return await dao.find_molecules_by_substructure(identifier)
