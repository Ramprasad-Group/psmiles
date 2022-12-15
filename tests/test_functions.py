from psmiles import PolymerSmiles as PS


def test_linear_copolymer():
    ps1 = PS("[*]CC[*]")
    ps2 = PS("[*]C=C[*]")
    lc = ps1.linear_copolymer(ps2)
    assert str(lc) == "[*]C=CC=CC=CCCCCCC[*]"


def test_block_copolymer():
    # Block polymer with 5A and 5B
    ps1 = PS("[*]CC[*]")
    ps2 = PS("[*]C=C[*]")
    lc = ps1.linear_copolymer(ps2, [0] * 5 + [1] * 5)
    assert str(lc) == "[*]C=CC=CC=CC=CC=CCCCCCCCCCC[*]"


def test_gradient_copolymer():
    # Gradient polymer
    ps1 = PS("[*]CC[*]")
    ps2 = PS("[*]C=C[*]")
    gradient_pattern = "AAAAAABAABBAABABBBAABBBBBB"
    lc = ps1.linear_copolymer(ps2, gradient_pattern)
    assert (
        str(lc)
        == "[*]C=CC=CC=CC=CC=CC=CCCCCC=CC=CC=CCCC=CCCCCC=CC=CCCCCC=CCCCCCCCCCCCC[*]"
    )


def test_random_copolymer():
    import random

    random.seed(10)
    ps1 = PS("[*]CC[*]")
    ps2 = PS("[*]CC([*])c1ccccc1")
    lc = ps1.random_copolymer(ps2, ratio=0.5, units=6)
    assert str(lc) == "[*]CC(CCCC(CCCCC(C[*])c1ccccc1)c1ccccc1)c1ccccc1"
