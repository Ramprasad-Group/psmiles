from __future__ import annotations

import importlib
import logging
import pprint
import re
import sys
from typing import Dict, List, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey
from sklearn.metrics.pairwise import cosine_similarity
from typing import Callable
from canonicalize_psmiles.canonicalize import canonicalize as ext_canonicalize


def copy_doc(copy_func: Callable) -> Callable:
    """Use Example: copy_doc(self.copy_func)(self.func) or used as deco"""

    def wrapper(func: Callable) -> Callable:
        func.__doc__ = copy_func.__doc__
        return func

    return wrapper


class PolymerSmiles:
    def __init__(self, psmiles: str, deactivate_warnings: bool = False):
        """PolymerSmiles ::  Fun with polymer SMILES strings.

        Note:
            Polymer SMILES strings must have two \* or [\*]
            Or, ladder polymers must have
            [e] <-> [t] and
            [d] <-> [g]

        Args:
            psmiles (str): PSMILES string, e.g. [*]CC[*]
            deactivate_warnings (bool, optional): Deactivate warnings. Defaults to False.
        """
        self.psmiles = psmiles
        self.ladder = False
        self.generator_circular = rdFingerprintGenerator.GetMorganGenerator()

        # convert * to [*]
        stars_no_bracket = re.findall(r"(?<!\[)\*(?!\])", self.psmiles)
        if len(stars_no_bracket) == 2:
            self.psmiles = self.psmiles.replace("*", "[*]")

        # linear homopolymer
        ct_stars = self.psmiles.count("[*]")

        # ladder polymers
        ladder_et = self.psmiles.count("[e]") + self.psmiles.count("[t]")
        ladder_dg = self.psmiles.count("[d]") + self.psmiles.count("[g]")

        if not deactivate_warnings:
            assert (
                ct_stars == 2 or ladder_dg == 2 or ladder_et == 2
            ), f"Smiles must have two [*], two *, [e] and [t], or [d] and [g] : {self.psmiles}"

        # Check if ladder PSMILES string
        if ladder_et == 2 or ladder_dg == 2:
            self.ladder = True

        # Check
        if not self.ladder:
            m = Chem.MolFromSmiles(self.psmiles, sanitize=False)
            if m is None:
                raise UserWarning(f"Invalid SMILES string: {self.psmiles}")
            else:
                error = Chem.SanitizeMol(m, catchErrors=True)
                if error:
                    raise UserWarning(
                        f"Invalid chemistry of {self.psmiles}. Issue with {error}"
                    )
        if not self.ladder and not deactivate_warnings:
            # Check double bonds
            self.check_double_bonds_at_connection()

        if self.ladder:
            logging.warning(
                "Ladder polymer detected. Only PG fingerprints are "
                "tested for ladder polymers."
            )

    def check_double_bonds_at_connection(self):
        """Check if SMILES has double bonds at both *
        (and not only at on *)
        """

        # get connection info
        info = self.get_connection_info()

        if info["neighbor"]["bond_type"][0] != info["neighbor"]["bond_type"][1]:
            raise UserWarning(
                f"The bond types of the SMILES string {self.psmiles} "
                f"at the connection points (*) is not the same."
                f"Bond types: {info['neighbor']['bond_type'][0]} - {info['neighbor']['bond_type'][1]}"
            )

    def get_connection_info(self, mol: Chem.RWMol = None, symbol="*") -> Dict:
        """Get connection information of stars and neighbors.

        Args:
            mol (Chem.RWMol, optional): _description_. Defaults to None.
            symbol (str, optional): _description_. Defaults to "*".

        Returns:
            Dict: Dictionary with some info
        """
        ret_dict = {}
        if mol is None:
            mol = self.mol

        stars_indices, stars_type, all_symbols, all_index = [], [], [], []
        for star_idx, atom in enumerate(mol.GetAtoms()):
            all_symbols.append(atom.GetSymbol())
            all_index.append(atom.GetIdx())
            if symbol in atom.GetSymbol():
                stars_indices.append(star_idx)
                stars_type.append(atom.GetSmarts())

        stars_bond = mol.GetBondBetweenAtoms(stars_indices[0], stars_indices[1])
        if stars_bond:
            stars_bond = stars_bond.GetBondType()

        ret_dict["symbols"] = all_symbols
        ret_dict["index"] = all_index

        ret_dict["star"] = {
            "index": stars_indices,
            "atom_type": stars_type,
            "bond_type": stars_bond,
        }

        # multiple neighbors are possible
        neighbor_indices = [
            [x.GetIdx() for x in mol.GetAtomWithIdx(stars_indices[0]).GetNeighbors()],
            [x.GetIdx() for x in mol.GetAtomWithIdx(stars_indices[1]).GetNeighbors()],
        ]

        neighbors_type = [
            [mol.GetAtomWithIdx(x).GetSmarts() for x in neighbor_indices[0]],
            [mol.GetAtomWithIdx(x).GetSmarts() for x in neighbor_indices[1]],
        ]

        # Bonds between stars and neighbors
        neighbor_bonds = [
            [
                mol.GetBondBetweenAtoms(stars_indices[0], x).GetBondType()
                for x in neighbor_indices[0]
            ],
            [
                mol.GetBondBetweenAtoms(stars_indices[1], x).GetBondType()
                for x in neighbor_indices[1]
            ],
        ]
        s_path = None
        if neighbor_indices[0][0] != neighbor_indices[1][0]:
            s_path = Chem.GetShortestPath(
                mol, neighbor_indices[0][0], neighbor_indices[1][0]
            )

        ret_dict["neighbor"] = {
            "index": neighbor_indices,
            "atom_type": neighbors_type,
            "bond_type": neighbor_bonds,
            "path": s_path,
        }

        # Stereo info
        stereo_info = []
        for b in mol.GetBonds():
            bond_type = b.GetStereo()
            if bond_type != Chem.rdchem.BondStereo.STEREONONE:
                idx = [b.GetBeginAtomIdx(), b.GetEndAtomIdx()]
                neigh_idx = b.GetStereoAtoms()
                stereo_info.append(
                    {
                        "bond_type": bond_type,
                        "atom_idx": idx,
                        "bond_idx": b.GetIdx(),
                        "neighbor_idx": list(neigh_idx),
                    }
                )

        ret_dict["stereo"] = stereo_info

        # Ring info
        ring_info = mol.GetRingInfo()
        ret_dict["atom_rings"] = ring_info.AtomRings()
        ret_dict["bond_rings"] = ring_info.BondRings()

        return ret_dict

    def __repr__(self) -> str:
        st = f"{self.psmiles}"
        return st

    def __str__(self) -> str:
        return self.psmiles

    def _repr_png_(self):
        """Display PSMILES PNG in Jupyter notebook."""
        print(f"SMILES: {self.__repr__()}")
        if not self.ladder:
            return self.mol._repr_png_()

    def replace_stars(self, _with: str) -> PolymerSmiles:
        """Replace stars with other characters

        Args:
            _with (str): Replacement characters

        Returns:
            PolymerSmiles: replaced polymer SMILES
        """
        return PolymerSmiles(
            self.psmiles.replace("[*]", _with), deactivate_warnings=True
        )

    @property
    def randomize(self) -> PolymerSmiles:
        """Return randomized PSMILES string

        Returns:
            PolymerSmiles: randomized PSMILES string
        """
        sm = Chem.MolToSmiles(
            Chem.MolFromSmiles(self.psmiles), doRandom=True, canonical=False
        )
        sm = sm.replace("*", "[*]")
        return PolymerSmiles(sm)

    def nb_display(self, mol):
        """Helper function to display polymer for debug"""
        print(f"SMILES: {Chem.MolToCXSmiles(mol)}")
        display(mol)

    @property
    def cyclic(self) -> PolymerSmiles:
        """Returns a cyclic SMILES string

        Returns:
            PolymerSmiles: cyclic SMILES string
        """
        logging.warning("Function is experimental. Please check results carefully.")

        mol = Chem.RWMol(Chem.MolFromSmiles(self.psmiles))

        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        atom_idx_star = [n for n, sym in enumerate(symbols) if sym == "*"]

        # Chose bondtyp ~
        bond_type = Chem.rdchem.BondType.UNSPECIFIED

        mol.AddBond(atom_idx_star[0], atom_idx_star[1], bond_type)
        sm = Chem.MolToSmiles(mol)
        sm = sm.replace("*", "[*]")

        return PolymerSmiles(sm)

    @property
    @copy_doc(ext_canonicalize)
    def canonicalize(self) -> PolymerSmiles:
        return PolymerSmiles(ext_canonicalize(self.psmiles))

    @property
    def inchi(self) -> str:
        """Returns the InChI string of the SMILES.

        [\*] is replaced with [At] in the output. Uses RDKit's MolToInchi

        Returns:
            str: InChI string
        """
        return MolToInchi(Chem.MolFromSmiles(self.psmiles.replace("[*]", "[At]")))

    @property
    def inchi_key(self) -> str:
        """Returns the InChI key of the SMILES.

        [\*] is replaced with [At]. Uses RDKit's MolToInchiKey

        Returns:
            str: InChI key
        """
        return MolToInchiKey(Chem.MolFromSmiles(self.psmiles.replace("[*]", "[At]")))

    @property
    def as_dict(self) -> dict:
        """Returns a dictionary of all transformations.

        Inchi/keys are computed from the canonicalized SMILES string

        Returns:
            dict: dictionary of SMILES transformations
        """
        d = dict(
            smiles=self.psmiles,
            dimer=self.dimer.__str__(),
            canonicalized=self.canonicalize.__str__(),
            inchi=self.canonicalize.inchi,
            inchi_key=self.canonicalize.inchi_key,
        )
        return d

    @property
    def dimer(self) -> PolymerSmiles:
        """Dimerize the SMILES string

        Returns:
            PolymerSmiles: dimerized SMILES string
        """
        # Make atom indices visable
        if logging.DEBUG >= logging.root.level:
            from rdkit.Chem.Draw import IPythonConsole

            IPythonConsole.drawOptions.addAtomIndices = True

        _m = Chem.MolToSmiles

        mol = self.mol
        info = self.get_connection_info(mol)
        logging.debug(f"(1) Get connection info \n {pprint.pformat(info)}")
        if logging.DEBUG >= logging.root.level:
            self.nb_display(mol)

        # combine two mols
        logging.debug(f"(2) Combine two mols")

        mol_combined = Chem.RWMol(Chem.CombineMols(mol, mol))
        if logging.DEBUG >= logging.root.level:
            self.nb_display(mol_combined)

        # Connect with single always
        bond_type = Chem.rdchem.BondType.SINGLE

        # Remove stars and add bonds between neighbors

        # Two connection possibilities
        connect = [
            info["star"]["index"][0],
            info["star"]["index"][1] + len(info["symbols"]),
        ]
        connect2 = [
            info["star"]["index"][0],
            info["star"]["index"][0] + len(info["symbols"]),
        ]

        logging.debug(
            f"(3) Connect star atoms {connect[0]} and {connect[1]} with {bond_type = }"
        )
        logging.debug(
            f"Alternative: connect star atoms {connect2[0]} and {connect2[1]} with {bond_type = }"
        )

        mol_combined.AddBond(connect[0], connect[1], order=bond_type)
        if logging.DEBUG >= logging.root.level:
            self.nb_display(mol_combined)

        # Remove patterns
        sm = Chem.MolToSmiles(mol_combined)
        patterns = {
            "**": "",  # normal bond
            "//": "/",  # if * have stereochemistry
            "\\\\": "\\",  # if * have stereochemistry
            "==": "=",  # if * has double bonds
            "##": "#",  # if * has triple bonds
            "\/": "\\",  # if \ or / at * (at the double bond)
        }
        logging.debug(f"(4) Remove {patterns} pattern {sm}")
        for p, r in patterns.items():
            sm = sm.replace(p, r)
            logging.debug(f"Replacing {p} with {r}: {sm}")

        if logging.DEBUG >= logging.root.level:
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(sm)))
            self.nb_display(mol)

        # Get dimer smiles
        # logging.debug(f"(5) Renumber atoms and get dimer smiles")
        # sm = Chem.MolToSmiles(mol_combined)

        return PolymerSmiles(sm)

    @property
    def mol(self) -> Chem.RWMol:
        """Returns a RDKit mol object.

        Note:
            In jupyter notebooks, this function draws the SMILES string

        Returns:
            Chem.MolFromSmiles: Mol object
        """
        return Chem.RWMol(Chem.MolFromSmiles(self.psmiles))

    @property
    def show_all(self) -> Chem.Draw.MolsToGridImage:
        """Draws all SMILES string variants and plot

        Note:
            Does only work in jupyter notebooks.

        Returns:
            Chem.Draw.MolsToGridImage: Drawing
        """

        mols = [
            self.mol,
            self.canonicalize.mol,
            self.dimer.mol,
        ]
        names = [
            "SMI=" + self.psmiles,
            "Can= " + self.canonicalize.psmiles,
            "Dim=" + self.dimer.psmiles,
        ]

        return Chem.Draw.MolsToGridImage(
            mols, molsPerRow=4, legends=names, subImgSize=(250, 200)
        )

    def fingerprint(self, fp="PG"):
        """Returns fingerprints of the PSMILES string.

        Note:
            PSMILES strings are canonicalized for the computation
            of the CI, mordred, and RDKit fingerprints.

        Args:
            fp (str, optional): Choose fingerprint from PG, CI, rdkit, or mordred. Defaults to 'PG'.

        Returns:
            _type_: Fingerprint vector
        """
        if fp == "PG":
            return self.fingerprint_pg
        elif fp == "CI":
            return self.fingerprint_circular
        elif fp == "mordred":
            return self.fingerprint_mordred
        elif fp == "rdkit":
            return self.fingerprint_rdkit
        elif fp == "polyBERT":
            return self.fingerprint_polyBERT
        else:
            raise UserWarning(f"Fingerprint {fp} unknown.")

    @property
    def fingerprint_polyBERT(self) -> np.ndarray:
        """Returns the PG fingerprint of the PSMILES string

        Returns:
            np.ndarray: polyBERT fingerprints
        """
        assert importlib.util.find_spec("sentence_transformers"), (
            "sentence-transformers python package is not installed. "
            "Please install with `poetry install --with polyBERT."
        )

        from sentence_transformers import SentenceTransformer

        polyBERT = SentenceTransformer('kuelumbus/polyBERT')

        return polyBERT.encode([self.canonicalize.psmiles])

    @property
    def fingerprint_pg(self) -> Dict[str, float]:
        """Returns the PG fingerprint of the PSMILES string

        Returns:
            Dict[str, float]: PG fingerprints
        """
        assert importlib.util.find_spec("pgfingerprinting"), (
            "pgfingerprinting python package is not installed. "
            "Please install with `poetry install --with pg. "
            "Package not available to the public."
        )

        from pgfingerprinting import fp as pgfp

        return pgfp.fingerprint_from_smiles(self.psmiles)

    @property
    def fingerprint_mordred(self) -> Dict[str, float]:
        """Return mordred fingerprints

        Returns:
            Dict[str, float]: mordred fingerprints
        """
        assert importlib.util.find_spec(
            "mordred"
        ), "Mordred is not installed. Please install with `poetry install --with mordred`"
        from mordred import Calculator, descriptors

        calc = Calculator(descriptors, ignore_3D=True)
        cpsmiles = self.canonicalize
        dim = calc.pandas([cpsmiles.dimer.replace_stars("[At]").mol], quiet=True, nproc=1)
        mon = calc.pandas([cpsmiles.replace_stars("[At]").mol], quiet=True, nproc=1)
        fps = dim.fill_missing().T - mon.fill_missing().T

        return fps[0].to_dict()

    @property
    def fingerprint_circular(self) -> np.ndarray:
        """Returns the circular (Morgen) count fingerprint

        Note that these fingerprint vary with the writing of
        the PSMILES string. Canonicalize before computing the fingerprint.

        Returns:
            numpy.ndarray: circular fingerprint
        """

        return self.generator_circular.GetCountFingerprintAsNumPy(self.canonicalize.mol).astype(int)

    @property
    def fingerprint_rdkit(self) -> np.ndarray:
        """Returns the RDKit count fingerprint

        Returns:
            numpy.ndarray: RDKit fingerprint
        """
        from rdkit.Chem import rdFingerprintGenerator

        fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator()
        fp_mono = fp_gen.GetCountFingerprintAsNumPy(self.canonicalize.mol).astype(int)
        return fp_mono

    def is_similar(self, other: PolymerSmiles) -> float:
        """Measures the cosine similarity between two SMILES stings.

        Args:
            other (PolymerSmiles): other SMILES string

        Returns:
            float: cosine similarity
        """

        fp1 = self.fingerprint
        fp2 = other.fingerprint
        intersection = [[fp1[x], fp2[x]] for x in fp1 if x in fp2]
        intersection = np.array(intersection).T
        return round(cosine_similarity(intersection)[0, 1], 5)

    def alternating_copolymer(self, other: Union[PolymerSmiles, str]):
        """Combines two SMILES strings to one (composition 0.5) alternating copolymer,
           essentially a homopolymer.

        TODO:
            Should have a switch to define how to combine two polymers

        Args:
            other (PolymerSmiles): Second polymer string

        Returns:
            PolymerSmiles: copolymer SMILES
        """
        logging.warning("Function is experimental. Please check results carefully.")

        if not isinstance(other, PolymerSmiles):
            other = PolymerSmiles(other)

        _m = Chem.MolToSmiles

        symbols1 = [a.GetSymbol() for a in self.mol.GetAtoms()]
        idx_star1 = [n for n, sym in enumerate(symbols1) if sym == "*"]

        symbols2 = [a.GetSymbol() for a in other.mol.GetAtoms()]
        idx_star2 = [n for n, sym in enumerate(symbols2) if sym == "*"]

        # combine two mols
        ed = Chem.RWMol(Chem.CombineMols(self.mol, other.mol))

        # Chose bondtyp ~
        bond_type = Chem.rdchem.BondType.UNSPECIFIED

        # Connect stars
        # Can be [0,0], [0,1], [1,0], [1,1]
        ed.AddBond(idx_star1[0], idx_star2[0] + len(symbols1), order=bond_type)

        # Get dimer smiles
        sm = _m(ed)

        # Remove stars connection *~*, and cis and trans around *~*
        patterns = [r"\\*~*\\", "/*~*/", "/*~*\\", "\\*~*/", "*~*"]
        for pat in patterns:
            sm = sm.replace(pat, "")

        return PolymerSmiles(sm)
