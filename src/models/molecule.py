from sqlalchemy import Column, Integer, String, Text
from sqlalchemy.orm import declarative_base

Base = declarative_base()


class Molecule(Base):
    __tablename__ = "molecules"
    id = Column(Integer, primary_key=True, index=True)
    identifier = Column(String, unique=True, index=True)
    smiles = Column(Text, nullable=False)
