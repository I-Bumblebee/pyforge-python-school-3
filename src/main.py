import pytest
from rdkit import Chem


def substructure_search(mols, mol):
    substructure = Chem.MolFromSmiles(mol)
    return list(filter(lambda smiles: Chem.MolFromSmiles(smiles).HasSubstructMatch(substructure), mols))


@pytest.mark.parametrize("mols, substructure, expected", [
    (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
     "c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
])
def test_substructure_search(mols, substructure, expected):
    assert substructure_search(mols, substructure) == expected
