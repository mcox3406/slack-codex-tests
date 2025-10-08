"""Utility to highlight rdchiral atom-map reassignment and a simple fix."""

from rdkit import Chem

from rdchiral import initialization
from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRun


reactant_smiles = '[CH2:1]([CH:2]=[CH:3][CH2:4][Cl:900])[CH:7]1[CH2:6][NH:5][CH:9]=[CH:8]1'
mech = '[C;H1;+0:4]-[C;+0:5]-[C;+0:3]-[C;+0:2]=[C;H1;+0:1]>>[C;H2;+0:1]-[C;+0:2]=[C;+0:3].[C;H0;+0:4]=[C;+0:5]'


def print_atom_maps(header: str, mol: Chem.Mol) -> None:
    """Print atom indices, symbols, and map numbers for a molecule."""

    print(header)
    for atom in mol.GetAtoms():
        print(atom.GetIdx(), atom.GetSymbol(), atom.GetAtomMapNum())


def run_reaction(label: str) -> None:
    """Execute the reaction and print the outcomes with keep_mapnums=True."""

    reactants = rdchiralReactants(reactant_smiles)
    print_atom_maps(
        f"\nAtom maps after rdchiralReactants initialization ({label}):",
        reactants.reactants,
    )
    rxn = rdchiralReaction(mech)
    products = rdchiralRun(rxn, reactants, keep_mapnums=True)
    print(f"\nProducts with keep_mapnums=True ({label}):")
    for prod in products:
        print(prod)


# Show the default behaviour, which overwrites the incoming map numbers
original = Chem.MolFromSmiles(reactant_smiles)
print_atom_maps('Original atom map numbers:', original)
run_reaction('vanilla rdchiral')


def initialize_reactants_from_smiles_preserve_maps(smiles: str) -> Chem.Mol:
    """Drop-in replacement that preserves existing atom maps from the input."""

    mol = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
    mol.UpdatePropertyCache(strict=False)
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            atom.SetAtomMapNum(atom.GetIdx() + 1)
    return mol


# Monkey patch rdchiral to use the preserving variant for subsequent runs
initialization.initialize_reactants_from_smiles = initialize_reactants_from_smiles_preserve_maps

run_reaction('patched to preserve input maps')
