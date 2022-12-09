from psmiles import PolymerSmiles as PS


def test_dimer_double_bonds_at_star():
    sms = [
        "N(=[*])c1ccc(NC=[*])cc1",
        "[*]=CNc1ccc(cc1)N=[*]",
        "[*]=Nc1ccc(NC=[*])cc1",
        "[*]=Nc1ccc(cc1)NC=[*]",
        "c1(ccc(NC=[*])cc1)N=[*]",
        "c1c(NC=[*])ccc(c1)N=[*]",
        "c1c(ccc(NC=[*])c1)N=[*]",
        "c1cc(N=[*])ccc1NC=[*]",
        "c1cc(N=[*])ccc1NC=[*]",
        "c1cc(ccc1N=[*])NC=[*]",
    ]

    dimerized = []
    for s in sms:
        dimerized.append(PS(s).canonicalize.dimer().psmiles)
    one_smiles = list(set(dimerized))

    assert one_smiles == ["[*]c1ccc(N=CNNC=Nc2ccc([*])cc2)cc1"]


def test_dimer_double_bonds_at_star_complicated():
    sms = [
        "C(=[*])(Nc1c(cc(cc1)-c1ccc(NC(=Nc2ccc(N=[*])cc2)C)c(C(=O)O)c1)C(=O)O)C",
        "C(Nc1c(cc(cc1)-c1cc(C(O)=O)c(NC(=[*])C)cc1)C(O)=O)(=Nc1ccc(cc1)N=[*])C",
        "CC(Nc1c(cc(cc1)-c1cc(C(O)=O)c(NC(C)=Nc2ccc(cc2)N=[*])cc1)C(=O)O)=[*]",
        "N(=[*])c1ccc(cc1)N=C(C)Nc1ccc(cc1C(O)=O)-c1ccc(c(C(=O)O)c1)NC(C)=[*]",
        "N(c1c(C(=O)O)cc(cc1)-c1ccc(c(C(=O)O)c1)NC(=[*])C)C(=Nc1ccc(N=[*])cc1)C",
        "c1(c(ccc(-c2cc(C(O)=O)c(NC(=[*])C)cc2)c1)NC(C)=Nc1ccc(N=[*])cc1)C(O)=O",
        "c1c(-c2cc(c(cc2)NC(=Nc2ccc(cc2)N=[*])C)C(O)=O)ccc(NC(=[*])C)c1C(=O)O",
        "c1cc(-c2ccc(NC(=[*])C)c(c2)C(O)=O)cc(C(=O)O)c1NC(C)=Nc1ccc(N=[*])cc1",
        "c1cc(-c2ccc(NC(C)=[*])c(c2)C(=O)O)cc(c1NC(=Nc1ccc(cc1)N=[*])C)C(=O)O",
        "c1cc(N=[*])ccc1N=C(Nc1c(cc(-c2ccc(c(C(=O)O)c2)NC(=[*])C)cc1)C(O)=O)C",
    ]

    dimerized = []
    for s in sms:
        dimerized.append(PS(s).canonicalize.dimer().psmiles)
    one_smiles = list(set(dimerized))

    assert one_smiles == [
        "[*]=C(C)Nc1ccc(-c2ccc(NC(C)=Nc3ccc(N=Nc4ccc(N=C(C)Nc5ccc(-c6ccc(NC(=[*])C)c(C(=O)O)c6)cc5C(=O)O)cc4)cc3)c(C(=O)O)c2)cc1C(=O)O"
    ]


def test_dimer_stereo_chemistry():
    sms = [
        "C(\\C(CCCOC(=O)/C=C/c1ccccc1)C[*])=C\\[*]",
        "c1cccc(/C=C/C(=O)OCCCC(C[*])/C=C\\[*])c1",
        "C(C(/C=C\\[*])C[*])CCOC(/C=C/c1ccccc1)=O",
        "C(\\c1ccccc1)=C/C(OCCCC(C[*])/C=C\\[*])=O",
        "C(OC(=O)/C=C/c1ccccc1)CCC(/C=C\\[*])C[*]",
        "C(CCOC(=O)/C=C/c1ccccc1)C(/C=C\\[*])C[*]",
        "c1(/C=C/C(OCCCC(/C=C\\[*])C[*])=O)ccccc1",
        "C(CCCOC(=O)/C=C/c1ccccc1)(C[*])/C=C\\[*]",
        "c1(/C=C/C(=O)OCCCC(C[*])/C=C\\[*])ccccc1",
        "C(CCCOC(=O)/C=C/c1ccccc1)(C[*])/C=C\\[*]",
    ]

    dimerized = []
    for s in sms:
        dimerized.append(PS(s).canonicalize.dimer().psmiles)
    one_smiles = list(set(dimerized))

    assert one_smiles == [
        "[*]C(CC=CC=CCC([*])CCCOC(=O)/C=C/c1ccccc1)CCCOC(=O)/C=C/c1ccccc1"
    ]


def test_dimer_stereo_chemistry_cis():
    sms = [
        "C([*])C/C=C\\[*]",
        "[*]/C=C\\CC[*]",
        "C(\\[*])=C\\CC[*]",
        "[*]CC/C=C\\[*]",
        "C(\\CC[*])=C\\[*]",
        "C(=C\\CC[*])\\[*]",
        "[*]/C=C\\CC[*]",
        "C(C/C=C\\[*])[*]",
        "C(\\[*])=C\\CC[*]",
        "[*]/C=C\\CC[*]",
    ]

    dimerized = []
    for s in sms:
        dimerized.append(PS(s).canonicalize.dimer().psmiles)
    one_smiles = list(set(dimerized))

    assert one_smiles == ["[*]=CCCC=CCCC=[*]"]


def test_dimer_stereo_chemistry_trans():
    sms = ["[*]/C=C(/[*])C(C)(C)", r"[*]\C=C(\[*])C(C)(C)"]
    dimerized = []
    for s in sms:
        dimerized.append(PS(s).canonicalize.dimer().psmiles)
    one_smiles = list(set(dimerized))

    assert one_smiles == ["[*]/C(=C/C=C(/[*])C(C)C)C(C)C"]
