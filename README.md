# PSMILES - Fun with P🙂 strings

PSMILES (Polymer SMILES or P🙂) strings are string representations of polymer structures. PSMILES are built upon the SMILES chemical language. The `PSMILES` Python package contains tools to manipulate and handle PSMILES strings.


## PSMILES string
A PSMILES string has two stars (`[*]` or `*`) symbols that indicate the two endpoints of the polymer repeat unit, and otherwise follow the daylight SMILES syntax defined at [OpenSmiles](http://opensmiles.org/opensmiles.html). For ladder polymers, the repeat unit is indicated by `[e]` -> `[t]` and `[d]` -> `[g]`. See [PSMILES guide](https://www.polymergenome.org/guide/) for more details.

Examples:

 Polyethylene | Polyethylene oxide | Polypropylene |
|-|-|-|
| `[*]CC[*]` | `[*]CCO[*]` | `[*]CC([*])C` | 

## Features, functions, and roadmap

- [x] Canonicalization of PSMILES
- [x] Dimerization of PSMILES
- [x] Fingerprints (numerical representation)
    - [x] Polymer Genome fingerprints (Ramprasad group internal only, not available to the public)
    - [x] Mordred fingerprints [https://github.com/mordred-descriptor/mordred](https://github.com/mordred-descriptor/mordred)
    - [x] Circular (Morgen) fingerprints as implemented in RDKit
    - [x] RDKit fingerprints as implemented in RDKit
- [x] Randomize PSMILES
- [x] Polymers similarity based on fingerprints
- [x] Copolymers from two PSMILES
- [x] Ladder polymers (only PG fingerprints so far)
- [ ] More

## Canonicalization of PSMILES

The raw PSMILES syntax is ambiguous and non-unique; i.e., the same polymer may be represented using many PSMILES string:

Polyethylene | Polyethylene oxide | Polypropylene |
|-|-|-|
| `[*]C[*]`   | `[*]CCO[*]` | `[*]CC([*])C` | 
| `[*]CC[*]`  | `[*]COC[*]` | `[*]CC(CC([*])C)C` | 
| `[*]CCC[*]` | `[*]OCC[*]` | `CC([*])C[*]` | 

The canonicalization routine of the `PSMILES` packages finds a canonicalized version of the SMILES string by

1. Finding the shortest representation of a PSMILES string 

`[*]CCOCCO[*]` ->  `[*]CCO[*]`

2. Making the PSMILES string cyclic

`[*]CCO[*]` -> `C1 CCO C1`

3. Applying the canonicalization routine as implemented in RDKit

`C1 CCO C1` -> `C1 COC C1`

4. Breaking the cyclic bond

`C1 COC C1` -> `[*]COC[*]`


## Install with poetry 

```bash
poetry add git+ssh://git@github.com/Ramprasad-Group/psmiles.git
```

With PG (not available to the public) and mordred fingerprints

```bash
poetry add git+ssh://git@github.com/Ramprasad-Group/psmiles.git --with pg --with mordred
```

## Install for development


1. Clone project
```sh
git clone git@github.com:Ramprasad-Group/psmiles.git
cd psmiles
poetry config virtualenvs.in-project true
poetry install -E pg -E mordred
```

## Usage

[`test.ipynb`](tests/test.ipynb) - shows the general usage of the package


