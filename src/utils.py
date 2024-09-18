from functools import wraps
from typing import Any, AsyncIterator, Callable, TypeVar

from rdkit import Chem

from models.molecule import Molecule

T = TypeVar('T')


def limit_async_generator(
    async_gen_func: Callable[..., AsyncIterator[T]],
) -> Callable[..., AsyncIterator[T]]:
    @wraps(async_gen_func)
    async def wrapper(*args: Any, **kwargs: Any) -> AsyncIterator[T]:
        limit = kwargs.get('limit', None)

        async for item in async_gen_func(*args, **kwargs):
            if limit is not None and limit <= 0:
                break
            yield item
            if limit is not None:
                limit -= 1

    return wrapper


def substructure_matches(
    molecule: Molecule, mol_to_match_to: Chem.Mol | None
) -> bool:
    """
    Checks if the given molecule matches the substructure of \
        the provided reference molecule.

    Parameters:
    - molecule (Molecule): An instance of the `Molecule` model. \
        The molecule to be checked.
    - mol_to_match_to (Chem.Mol | None): An instance of `Chem.Mol` \
        representing the reference molecule. \
        If `None`, the function will return `True`.

    Returns:
    - bool: `True` if either:
        - `mol_to_match_to` is `None` (indicating no reference molecule to \
            match against), or
        - The `molecule` has a substructure match with `mol_to_match_to`.
      Returns `False` if:
        - `mol_to_match_to` exists and the `molecule` does not match \
            the substructure of `mol_to_match_to`.
    """
    if mol_to_match_to is None:
        return True

    mol_to_check = Chem.MolFromSmiles(molecule.smiles)
    return mol_to_check.HasSubstructMatch(mol_to_match_to)
