from rdkit import Chem


def substructure_search(mols, mol):
    substructure = Chem.MolFromSmiles(mol)
    return [
        m for m in mols
        if (mol_to_check := Chem.MolFromSmiles(m.smiles)) is not None and
        mol_to_check.HasSubstructMatch(substructure)
    ]
