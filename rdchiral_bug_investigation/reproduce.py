from rdkit import Chem
from rdchiral.main import rdchiralRun, rdchiralReaction, rdchiralReactants

reactant_smiles = '[CH2:1]([CH:2]=[CH:3][CH2:4][Cl:900])[CH:7]1[CH2:6][NH:5][CH:9]=[CH:8]1'
mech = '[C;H1;+0:4]-[C;+0:5]-[C;+0:3]-[C;+0:2]=[C;H1;+0:1]>>[C;H2;+0:1]-[C;+0:2]=[C;+0:3].[C;H0;+0:4]=[C;+0:5]'

print('Original atom map numbers:')
original = Chem.MolFromSmiles(reactant_smiles)
for atom in original.GetAtoms():
    print(atom.GetIdx(), atom.GetSymbol(), atom.GetAtomMapNum())

reactants = rdchiralReactants(reactant_smiles)
print('\nAtom map numbers after rdchiralReactants initialization:')
for atom in reactants.reactants.GetAtoms():
    print(atom.GetIdx(), atom.GetSymbol(), atom.GetAtomMapNum())

rxn = rdchiralReaction(mech)
products = rdchiralRun(rxn, reactants, keep_mapnums=True)
print('\nProducts with keep_mapnums=True:')
for prod in products:
    print(prod)
