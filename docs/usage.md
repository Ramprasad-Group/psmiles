# Usage

Simply create an object of the class `PolymerSmiles` using your PSMILES string. The object has functions to compute properties or manipulate the PSMILES string. For example, canonicalize a PSMILES string with

```python
from psmiles import PolymerSmiles as PS

ps = PS("C(c1ccccc1)(C[*])[*]")
ps.canonicalize
```

Gives

``` py
PSMILES: [*]CC([*])c1ccccc1
```

If you work in a Jupyter notebook, psmiles returns you the chemical drawing too.


Get the polyBERT fingerprint with

```python
from psmiles import PolymerSmiles as PS

ps = PS("C(c1ccccc1)(C[*])[*]")
ps.fingerprint("polyBERT")
```

See the [test_book.ipynb](https://github.com/Ramprasad-Group/psmiles/blob/main/tests/test_book.ipynb) at GitHub or directly open it in [Colab](https://colab.research.google.com/github/Ramprasad-Group/psmiles/blob/main/tests/test_book.ipynb) for more examples.
