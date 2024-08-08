#!/usr/bin/env python

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy_utils import create_database, database_exists

from models import Base, Molecule

DATABASE_URL = "sqlite:///./src/db_data/smile.sqlite"

# Connect to database and create neccessary tables/database if not exists
engine = create_engine(DATABASE_URL)

allready_migrated = database_exists(engine.url)
if not allready_migrated:
    create_database(engine.url)

Base.metadata.create_all(engine)

if not allready_migrated:
    # Configure sessionmake
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

    # Start transaction and seed the database
    with SessionLocal() as db:
        molecules = [
            Molecule(identifier="mol1", smiles="C1=CC=CC=C1"),
            Molecule(identifier="mol2", smiles="C2H6"),
        ]
        db.add_all(molecules)
        db.commit()
