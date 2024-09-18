import logging

import pytest
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select

from configs.setup_logger import setup_logger
from src.models.molecule import Molecule

logger = logging.getLogger(__name__)
setup_logger()


@pytest.mark.asyncio
async def test_create_and_verify_molecule_in_db(
    async_client: AsyncClient,
    test_session: AsyncSession,
):
    molecule_data = {
        "identifier": "test123",
        "smiles": "C1=CC=CC=C1",
    }

    response = await async_client.post("/molecules/", json=molecule_data)

    assert response.status_code == 200

    created_molecule = (
        await test_session.execute(
            select(Molecule).filter_by(identifier=molecule_data["identifier"])
        )
    ).scalar_one()

    assert created_molecule is not None
    assert created_molecule.identifier == molecule_data["identifier"]
    assert created_molecule.smiles == molecule_data["smiles"]


@pytest.mark.asyncio
async def test_get_molecule(
    async_client: AsyncClient,
    test_session: AsyncSession,
):
    molecule = Molecule(
        identifier="mol_test",
        smiles="CCO",
    )
    test_session.add(molecule)

    response = await async_client.get(f"/molecules/{molecule.identifier}")

    assert response.status_code == 200

    data = response.json()
    assert data["identifier"] == molecule.identifier
    assert data["smiles"] == molecule.smiles


@pytest.mark.asyncio
async def test_update_molecule(
    async_client: AsyncClient,
    test_session: AsyncSession,
):
    molecule = Molecule(
        identifier="mol_update",
        smiles="C2H6",
    )
    test_session.add(molecule)

    new_data = {"identifier": "mol_update", "smiles": "CH4"}
    response = await async_client.put(
        f"/molecules/{molecule.identifier}",
        json=new_data,
    )

    assert response.status_code == 200

    updated_molecule = (
        await test_session.execute(
            select(Molecule).filter_by(identifier="mol_update")
        )
    ).scalar_one()

    assert updated_molecule is not None
    assert updated_molecule.smiles == new_data["smiles"]


@pytest.mark.asyncio
async def test_delete_molecule(
    async_client: AsyncClient,
    test_session: AsyncSession,
):
    molecule = Molecule(
        identifier="mol_delete",
        smiles="CO2",
    )
    test_session.add(molecule)

    response = await async_client.delete(f"/molecules/{molecule.identifier}")

    assert response.status_code == 200

    result = (
        await test_session.execute(
            select(Molecule).filter_by(identifier="mol_delete")
        )
    ).scalar_one_or_none()

    assert result is None
