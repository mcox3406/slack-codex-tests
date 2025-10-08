# rdchiral Atom Mapping Investigation

## Reproduction
- `rdchiral_bug_investigation/reproduce.py` instantiates the reaction from the Slack example and runs it with `keep_mapnums=True`.
- Running the script shows that the chlorine atom (originally map 900) is remapped to 5 and other atoms are renumbered sequentially despite `keep_mapnums=True`.

```
$ python reproduce.py
Original atom map numbers:
...
4 Cl 900
...
Atom maps after rdchiralReactants initialization (vanilla rdchiral):
...
4 Cl 5
5 N 6
...
Products with keep_mapnums=True (vanilla rdchiral):
[CH2:1]=[C:2]=[CH:3][CH2:4][Cl:5].[CH:6]1=[CH:10][CH2:9][NH:8][CH2:7]1
[CH2:1]=[CH:2][CH2:3][CH2:4][Cl:5].[CH:6]1=[C:10]=[CH:9][NH:8][CH2:7]1

Atom maps after rdchiralReactants initialization (patched to preserve input maps):
...
4 Cl 900
5 N 5
...
Products with keep_mapnums=True (patched to preserve input maps):
[CH2:1]=[C:2]=[CH:3][CH2:4][Cl:900].[NH:5]1[CH2:6][CH:7]=[CH:8][CH2:9]1
[CH2:1]=[CH:2][CH2:3][CH2:4][Cl:900].[NH:5]1[CH2:6][CH:7]=[C:8]=[CH:9]1
```

## Root Cause
The renumbering happens during rdchiral's reactant initialization. The current implementation overwrites every atom map number instead of preserving the values supplied in the input SMILES:

```
# rdchiral/initialization.py
[a.SetAtomMapNum(i+1) for (i, a) in enumerate(reactants.GetAtoms())]
```

Because the original mapping is discarded before `rdchiralRun` executes, subsequent steps cannot recover the initial identifiers, leading to the incorrect output observed above.

## Suggested Fix
Adjust `initialize_reactants_from_smiles` so that it only assigns sequential atom map numbers to atoms that are missing a map, leaving any pre-assigned identifiers untouched. The patched block in `reproduce.py` demonstrates that preserving the original atom map numbers keeps the chlorine at 900 and the nitrogen at 5 while still providing map numbers for previously unmapped atoms.

## Verification
Running the script after the monkey patch shows that the chlorinated carbon retains map 900 and the ring nitrogen retains map 5 in both products, matching the user-supplied identifiers while still producing duplicate product structures for completeness.
