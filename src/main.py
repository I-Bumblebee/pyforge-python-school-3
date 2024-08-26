from os import getenv
from typing import Optional

from fastapi import Depends, FastAPI, HTTPException, Query
from pydantic import BaseModel
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from models import Molecule
from utils import substructure_search
from config import settings


engine = create_engine(settings.database_url)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

app = FastAPI()


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


class MoleculeCreate(BaseModel):
    identifier: str
    smiles: str


@app.post("/molecules/")
def store_molecule(molecule: MoleculeCreate, db: Session = Depends(get_db)):
    db_molecule = Molecule(
        identifier=molecule.identifier,
        smiles=molecule.smiles)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    return db_molecule


@app.get("/molecules/{identifier}")
def view_molecule(identifier: str, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(
        Molecule.identifier == identifier).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule


@app.put("/molecules/{identifier}")
def update_molecule(
    identifier: str, molecule: MoleculeCreate, db: Session = Depends(get_db)
):
    db_molecule = db.query(Molecule).filter(
        Molecule.identifier == identifier).first()
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    db_molecule.smiles = molecule.smiles
    db.commit()
    db.refresh(db_molecule)
    return db_molecule


@app.delete("/molecules/{identifier}")
def destroy_molecule(identifier: str, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(
        Molecule.identifier == identifier).first()
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(db_molecule)
    db.commit()
    return {"detail": "Molecule deleted"}


@app.get("/molecules/")
def index_molecules(
        identifier: Optional[str] = Query(
            None,
            description="Identifier of a molecule whose \
                substructure will be used to find\
              and match other molecules with similar substructures."),
        db: Session = Depends(get_db),
):
    if identifier:
        mol_to_match = (db.query(Molecule).filter(
            Molecule.identifier == identifier).first())
        if not mol_to_match:
            raise HTTPException(
                status_code=404,
                detail="Molecule with given identifier not found")

        all_molecules = db.query(Molecule).all()
        filtered_molecules = substructure_search(
            all_molecules, mol_to_match.smiles)
        return filtered_molecules

    molecules = db.query(Molecule).all()
    return molecules


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}
