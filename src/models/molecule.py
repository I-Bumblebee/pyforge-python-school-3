from sqlalchemy import Column, Integer, String
from sqlalchemy.orm import declarative_base

Base = declarative_base()


class Molecule(Base):
    __tablename__ = "molecules"
    id = Column(Integer, primary_key=True)
    identifier = Column(String, unique=True)
    smiles = Column(String)
