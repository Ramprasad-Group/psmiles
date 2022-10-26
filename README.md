# PSMILES - Fun with P🙂 strings

The `PSMILES` Python package contains tools to work with polymer SMILES (PSMILES) strings.

## PSMILES strings

PSMILES (P🙂) strings are string representations of polymer structures (e.g., `[*]CC[*]` for polyethylene). PSMILES are built upon the SMILES chemical language. A PSMILES string has two stars (`[*]` or `*`) symbols that indicate the two endpoints of the polymer repeat unit and otherwise follows the daylight SMILES syntax defined at [OpenSmiles](http://opensmiles.org/opensmiles.html). See [PSMILES guide](https://www.polymergenome.org/guide/) for more details.

Examples:

 Polyethylene | Polyethylene oxide | Polypropylene |
|-|-|-|
| `[*]CC[*]` | `[*]CCO[*]` | `[*]CC([*])C` | 
| ![](docs/PE.png) | ![](docs/PEO.png) | ![](docs/PP.png) | 


## Features, functions, and roadmap

- [x] Canonicalize PSMILES strings (via https://github.com/Ramprasad-Group/canonicalize_psmiles)
- [x] Polymer fingerprints (descriptors or features)
    - [x] polyBERT fingerprints (see [polyBERT](https://github.com/Ramprasad-Group/polyBERT) and [arXiv](https://arxiv.org/abs/2209.14803)) 
    - [x] Polymer Genome fingerprints (Ramprasad group internal only, not available to the public)
    - [x] [Mordred](https://github.com/mordred-descriptor/mordred) fingerprints
    - [x] Circular (Morgen) fingerprints as implemented in RDKit
    - [x] RDKit fingerprints as implemented in RDKit
- [x] Dimerize PSMILES strings
- [x] Randomize PSMILES strings
- [x] Compute polymer similarity based on the fingerprints
- [x] Create alternating copolymers from two PSMILES strings


## Install with poetry 

```bash
poetry add git+https://github.com/Ramprasad-Group/psmiles.git
```

With polyBERT and mordred fingerprints

```bash
poetry add git+https://github.com/Ramprasad-Group/psmiles.git -E polyBERT -E mordred
```


## Install for development


1. Clone project
```sh
git clone https://github.com/Ramprasad-Group/psmiles.git
cd psmiles
poetry config virtualenvs.in-project true
poetry install --with polyBERT,mordred
```

## Usage

[`test_book.ipynb`](tests/test_book.ipynb) - shows the general usage of the package


