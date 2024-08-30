from typing import List, Optional

from fastapi import Depends, HTTPException
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select

from configs.database import get_db_session
from src.models.molecule import Molecule
from src.utils import substructure_search


class MoleculeDAO:
    def __init__(self, session: AsyncSession = Depends(get_db_session)):
        self.session = session

    async def create_molecule(self, identifier: str, smiles: str) -> Molecule:
        new_molecule = Molecule(identifier=identifier, smiles=smiles)
        self.session.add(new_molecule)
        return new_molecule

    async def get_molecule_by_identifier(self, identifier: str) -> Molecule:
        result = await self.session.execute(select(Molecule).filter_by(identifier=identifier))
        molecule = result.scalars().first()
        if not molecule:
            raise HTTPException(status_code=404, detail="Molecule not found")
        return molecule

    async def update_molecule(self, identifier: str, smiles: str) -> Molecule:
        result = await self.session.execute(select(Molecule).filter_by(identifier=identifier))
        molecule = result.scalars().first()
        if not molecule:
            raise HTTPException(status_code=404, detail="Molecule not found")

        molecule.smiles = smiles
        return molecule

    async def delete_molecule(self, identifier: str) -> None:
        result = await self.session.execute(select(Molecule).filter_by(identifier=identifier))
        molecule = result.scalars().first()
        if not molecule:
            raise HTTPException(status_code=404, detail="Molecule not found")

        await self.session.delete(molecule)

    async def find_molecules_by_substructure(
            self, identifier: Optional[str] = None) -> List[Molecule]:
        if identifier:
            result = await self.session.execute(select(Molecule).filter_by(identifier=identifier))
            mol_to_match = result.scalars().first()

            if not mol_to_match:
                raise HTTPException(
                    status_code=404,
                    detail="Molecule with given identifier not found")

            all_molecules_result = await self.session.execute(select(Molecule))
            all_molecules = all_molecules_result.scalars().all()

            filtered_molecules = substructure_search(
                all_molecules, mol_to_match.smiles)
            return filtered_molecules

        all_molecules_result = await self.session.execute(select(Molecule))
        return all_molecules_result.scalars().all()
