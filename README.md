# PSMILES - Fun with P🙂 strings

PSMILES (Polymer SMILES or P🙂) strings are string representations of polymer structures (e.g., `[*]CC[*]` for polyethylene). PSMILES are built upon the SMILES chemical language. The `PSMILES` Python package contains tools to manipulate and handle PSMILES strings.


## PSMILES string
A PSMILES string has two stars (`[*]` or `*`) symbols that indicate the two endpoints of the polymer repeat unit and otherwise follow the daylight SMILES syntax defined at [OpenSmiles](http://opensmiles.org/opensmiles.html). For ladder polymers, the repeat unit is indicated by `[e]` -> `[t]` and `[d]` -> `[g]`. See [PSMILES guide](https://www.polymergenome.org/guide/) for more details.

Examples:

 Polyethylene | Polyethylene oxide | Polypropylene |
|-|-|-|
| `[*]CC[*]` | `[*]CCO[*]` | `[*]CC([*])C` | 

## Features, functions, and roadmap

- [x] Canonicalization of PSMILES (via https://github.com/Ramprasad-Group/canonicalize_psmiles)
- [x] Dimerization of PSMILES
- [x] Fingerprints (numerical representation)
    - [x] polyBERT fingerprints (see arXiv) 
    - [x] Polymer Genome fingerprints (Ramprasad group internal only, not available to the public)
    - [x] Mordred fingerprints [https://github.com/mordred-descriptor/mordred](https://github.com/mordred-descriptor/mordred)
    - [x] Circular (Morgen) fingerprints as implemented in RDKit
    - [x] RDKit fingerprints as implemented in RDKit
- [x] Fingerprints for ladder polymers only for PG fingerprints
- [x] Randomize PSMILES
- [x] Polymer similarity based on fingerprints
- [x] Create alternating copolymers from two PSMILES


## Install with poetry 

```bash
poetry add git+ssh://git@github.com/Ramprasad-Group/psmiles.git
```

With polyBERT fingerprints

```bash
poetry add git+ssh://git@github.com/Ramprasad-Group/psmiles.git --with polyBERT
```

With mordred fingerprints

```bash
poetry add git+ssh://git@github.com/Ramprasad-Group/psmiles.git --with mordred
```

## Install for development


1. Clone project
```sh
git clone git@github.com:Ramprasad-Group/psmiles.git
cd psmiles
poetry config virtualenvs.in-project true
poetry install
```

## Usage

[`test.ipynb`](tests/test.ipynb) - shows the general usage of the package


