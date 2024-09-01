import pytest

from src.utils import substructure_search


@pytest.mark.parametrize(
    "mols, substructure, expected", [
        (
            ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
            "c1ccccc1",
            ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"],
        ),
        (["CCO", "CC(=O)O", "CC(=O)O"], "c1ccccc1", []),
        ([], "c1ccccc1", []),
        (["CCO", "c1ccccc1"], "", []),
    ]
)
def test_substructure_search(mols, substructure, expected):
    mols = [type('Mol', (object, ), {'smiles': smiles}) for smiles in mols]
    result = substructure_search(mols, substructure)
    assert [m.smiles for m in result] == expected
