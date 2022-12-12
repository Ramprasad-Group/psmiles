# PSMILES - Fun with P🙂s strings

The `psmiles` Python package provides tools and functions for polymer SMILES (PSMILES or P🙂s) strings. PSMILES strings are a chemical language to represent polymers.


[Install](install.md){ .md-button .md-button--primary }
[Usage](usage.md){ .md-button .md-button--primary }
[Run on Colab](https://colab.research.google.com/github/Ramprasad-Group/psmiles/blob/main/tests/test_book.ipynb){ .md-button .md-button--primary }

## What is a PSMILES string?

PSMILES strings are string representations of polymer chemical structures. PSMILES strings are very useful for data-driven polymer discovery, design or prediction task.

A PSMILES string follows the daylight SMILES syntax defined at [OpenSmiles](http://opensmiles.org/opensmiles.html), but has two stars (`[*]` or `*`) that indicate the two endpoints of the polymer repeat unit. See [PSMILES guide](https://www.polymergenome.org/guide) for more details.

Example P🙂s:

 Polyethylene | Polyethylene oxide | Polypropylene |
|-|-|-|
| `[*]CC[*]` | `[*]CCO[*]` | `[*]CC([*])C` |
| ![](PE.png) | ![](PEO.png) | ![](PP.png) |

??? tip

    Create these figures using `psmiles`
    ``` py linenums="1"
    from psmiles import PolymerSmiles as PS

    psmiles_strings = ['[*]CC[*]', '[*]CCO[*]', '[*]CC([*])C']
    [PS(ps).savefig() for ps in psmiles_strings]
    ```

## Features

- :fontawesome-solid-fingerprint: __Polymer__ __Fingerprints__ Numerical representations of polymers that measure polymer similarity. They can be used for any polymer informatics task that requires numerical representations of polymers such as property predictions, polymer structure predictions (design tasks), ML-based synthesis assistants, etc. `psmiles` offers polyBERT, Circular (Morgen), Mordred, and RDKit fingerprints.
- :simple-canonical: __Canonicalize__ __PSMILES__ Find a unique representation of the PSMILES string. Useful for many informatics tasks.
- :arrow_double_up: __Dimerize__ __PSMILES__ Get the dimerized PSMILES string
- :material-more: __More__ Radomize, compute polymer similarity, alternating copolymers, save chemical drawings
