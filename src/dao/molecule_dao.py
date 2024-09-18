import logging.config
from typing import AsyncIterator, Optional

from fastapi import Depends, HTTPException
from rdkit import Chem
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select

from configs.database import get_db_session
from configs.setup_logger import setup_logger
from src.models.molecule import Molecule
from src.utils import limit_async_generator, substructure_matches

logger = logging.getLogger(__name__)
setup_logger()


class MoleculeDAO:
    def __init__(self, session: AsyncSession = Depends(get_db_session)):
        self.session = session

    async def create_molecule(self, identifier: str, smiles: str) -> Molecule:
        new_molecule = Molecule(identifier=identifier, smiles=smiles)
        self.session.add(new_molecule)
        return new_molecule

    async def get_molecule_by_identifier(self, identifier: str) -> Molecule:
        query = select(Molecule).filter_by(identifier=identifier)
        molecule = (await self.session.execute(query)).scalar_one_or_none()

        if not molecule:
            raise HTTPException(status_code=404, detail="Molecule not found")
        return molecule

    async def update_molecule(self, identifier: str, smiles: str) -> Molecule:
        molecule = await self.get_molecule_by_identifier(identifier)
        molecule.smiles = smiles
        return molecule

    async def delete_molecule(self, identifier: str) -> None:
        molecule = await self.get_molecule_by_identifier(identifier)
        await self.session.delete(molecule)

    @limit_async_generator
    async def list_molecules(
        self,
        limit: Optional[int] = None,
        identifier: Optional[str] = None,
    ) -> AsyncIterator[Molecule]:
        query = select(Molecule).execution_options(stream_results=True)

        mol_to_match_to = None
        if identifier:
            molecule_to_match_to = await self.get_molecule_by_identifier(
                identifier
            )
            mol_to_match_to = Chem.MolFromSmiles(molecule_to_match_to.smiles)

        async for molecule in await self.session.stream_scalars(query):
            if identifier and not substructure_matches(
                molecule,
                mol_to_match_to,
            ):
                continue
            yield molecule
