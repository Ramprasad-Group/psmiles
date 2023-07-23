import numpy as np
from psmiles.psmiles import PolymerSmiles
from rdkit import Chem
from importlib import util
from rdkit.Chem import rdFingerprintGenerator


class Smiles(PolymerSmiles):
    def __init__(self, smiles: str, deactivate_warnings: bool = True):
        self.smiles = smiles
        psmiles = self.smiles
        super().__init__(psmiles, deactivate_warnings)

    def can_molecule(self):
        mol = Chem.MolFromSmiles(self.psmiles)
        return Chem.MolToSmiles(mol)

    @property
    def fingerprint_polyBERT(self) -> np.ndarray:
        """Compute the polyBERT fingerprint

        Note:
            Calling this will pull polyBERT from the hugging face hub.

        Returns:
            np.ndarray: polyBERT fingerprints
        """
        assert util.find_spec("sentence_transformers"), (
            "PolyBERT fingerprints require the 'sentence-transformers' Python package."
            " Please install with "
            "`pip install 'psmiles[polyBERT]@git+https://github.com/"
            "Ramprasad-Group/psmiles.git'` "
            "Or "
            "`poetry add git+https://github.com/"
            "Ramprasad-Group/psmiles.git -E polyBERT` "
        )

        can_smiles = self.can_molecule()

        from sentence_transformers import SentenceTransformer

        polyBERT = SentenceTransformer("kuelumbus/polyBERT")

        return polyBERT.encode([can_smiles], show_progress_bar=False)[0]

    @property
    def fingerprint_circular(self) -> np.ndarray:
        """Compute the circular (Morgen) count fingerprint
        
        Returns:
            numpy.ndarray: circular fingerprint
        """

        fp_gen = rdFingerprintGenerator.GetMorganGenerator()
        print(Chem.MolFromSmiles(self.smiles))
        return fp_gen.GetCountFingerprintAsNumPy(
            Chem.MolFromSmiles(self.smiles)
        ).astype(int)
