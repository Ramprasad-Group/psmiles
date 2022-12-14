# Usage

Simply create an object of the class `PolymerSmiles` using your PSMILES string. The object has functions to compute properties or manipulate the PSMILES string. For example, canonicalize a PSMILES string with

```py
>>> from psmiles import PolymerSmiles as PS
>>> ps = PS("C(c1ccccc1)(C[*])[*]")
>>> ps.canonicalize
[*]CC([*])c1ccccc1
```

If you work in a Jupyter notebook, this will also show the chemical drawing.


Get the polyBERT fingerprint with

```py
>>> from psmiles import PolymerSmiles as PS
>>> ps = PS("C(c1ccccc1)(C[*])[*]")
>>> ps.fingerprint("polyBERT")
[fingerprint]
```

Get the two dimers of the PSMILES string

```py
>>> from psmiles import PolymerSmiles as PS
>>> ps = PS("C(c1ccccc1)(C[*])[*]")
>>> ps.dimer(0)
[*]C(CCC([*])c1ccccc1)c1ccccc1
>>> ps.dimer(1)
[*]CC(CC([*])c1ccccc1)c1ccccc1
```

Create an alternating copolymer

```py
>>> from psmiles import PolymerSmiles as PS
>>> ps1 = PS('[*]CC[*]')
>>> ps2 = PS('[*]CCO[*]')
>>> ps1.alternating_copolymer(ps2, [0,0])
[*]CCCCO[*]
```

See the [test_book.ipynb](https://github.com/Ramprasad-Group/psmiles/blob/main/tests/test_book.ipynb) at GitHub or directly open it in [Colab](https://colab.research.google.com/github/Ramprasad-Group/psmiles/blob/main/tests/test_book.ipynb) for more examples.
