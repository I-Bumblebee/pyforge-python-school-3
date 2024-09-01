import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from src.main import app, get_db
from src.models.molecule import Base, Molecule
from configs.settings import settings


@pytest.fixture(scope='module')
def test_db():
    engine = create_engine(
        settings.database_url, connect_args={'check_same_thread': False}
    )
    Base.metadata.create_all(engine)
    yield engine
    Base.metadata.drop_all(engine)


@pytest.fixture(scope='function')
def db_session(test_db):
    connection = test_db.connect()
    transaction = connection.begin()
    session = sessionmaker(bind=connection)()
    yield session
    session.close()
    transaction.rollback()
    connection.close()


@pytest.fixture
def client(db_session):
    app.dependency_overrides[get_db] = lambda: db_session
    return TestClient(app)


@pytest.fixture
def mock_molecule_data(db_session: Session):
    molecules = [
        Molecule(identifier="test_id_1", smiles="C1=CC=CC=C1"),
        Molecule(identifier="test_id_2", smiles="CCO"),
        Molecule(identifier="test_id_3", smiles="CC(=O)Oc1ccccc1C(=O)O"),
        Molecule(identifier="test_id_4", smiles="c1ccccc1"),
    ]
    db_session.add_all(molecules)
    db_session.commit()

    return molecules


def test_get_all_molecules(client, mock_molecule_data):
    response = client.get("/molecules")
    assert response.status_code == 200
    assert len(response.json()) == len(mock_molecule_data)


# This endpoint first gets moelcule with provided identifier
# then executes substructure search on every other molcules and returns result
# result allways containes at least one molecule because molecule is
# substructure of itself
def test_index_molecules_with_identifier(client, mock_molecule_data):
    response = client.get("/molecules?identifier=test_id_1")
    assert response.status_code == 200
    assert len(response.json()) > 0


def test_index_molecules_with_nonexistent_identifier(client):
    response = client.get("/molecules?identifier=nonexistent_id")
    assert response.status_code == 404
    assert response.json() == {
        "detail": "Molecule with given identifier not found"
    }
