"""Unit tests for PointGroup factory functions.

These tests verify that every factory function returns the correct Schönflies
name and group order.  No pynauty or geometry is required — these are pure
data-layer tests on the PointGroup dataclass.
"""

import pytest

from polygraph.algorithms.symmetry.point_groups import (
    C_I,
    C_S,
    PointGroup,
    cyclic,
    cyclic_horizontal,
    cyclic_pure,
    dihedral,
    dihedral_chiral,
    icosahedral,
    improper,
    octahedral,
    pyritohedral,
    tetrahedral,
)

# ---------------------------------------------------------------------------
# PointGroup dataclass behaviour
# ---------------------------------------------------------------------------


def test_point_group_equality():
    assert PointGroup("T_d", 24) == PointGroup("T_d", 24)


def test_point_group_inequality_by_name():
    assert PointGroup("T_d", 24) != PointGroup("T_h", 24)


def test_point_group_inequality_by_order():
    assert PointGroup("T", 12) != PointGroup("T", 24)


def test_point_group_hashable():
    s = {PointGroup("O_h", 48), PointGroup("O_h", 48)}
    assert len(s) == 1


def test_point_group_repr():
    pg = PointGroup("I_h", 120)
    assert repr(pg) == "PointGroup('I_h', order=120)"


# ---------------------------------------------------------------------------
# Polyhedral (exceptional) group factories — all three tetrahedral groups
# ---------------------------------------------------------------------------


def test_tetrahedral_full():
    pg = tetrahedral(full=True)
    assert pg.name == "T_d"
    assert pg.order == 24


def test_tetrahedral_chiral():
    pg = tetrahedral(full=False)
    assert pg.name == "T"
    assert pg.order == 12


def test_pyritohedral():
    pg = pyritohedral()
    assert pg.name == "T_h"
    assert pg.order == 24


def test_three_tetrahedral_groups_are_distinct():
    """T, T_d, T_h are distinct even though T_d and T_h share order 24."""
    assert tetrahedral(full=False) != tetrahedral(full=True)
    assert tetrahedral(full=True) != pyritohedral()
    assert tetrahedral(full=False) != pyritohedral()


def test_octahedral_full():
    pg = octahedral(full=True)
    assert pg.name == "O_h"
    assert pg.order == 48


def test_octahedral_chiral():
    pg = octahedral(full=False)
    assert pg.name == "O"
    assert pg.order == 24


def test_icosahedral_full():
    pg = icosahedral(full=True)
    assert pg.name == "I_h"
    assert pg.order == 120


def test_icosahedral_chiral():
    pg = icosahedral(full=False)
    assert pg.name == "I"
    assert pg.order == 60


# ---------------------------------------------------------------------------
# Dihedral group factories
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_dihedral_h(n):
    pg = dihedral(n, "h")
    assert pg.name == f"D_{n}h"
    assert pg.order == 4 * n


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_dihedral_d(n):
    pg = dihedral(n, "d")
    assert pg.name == f"D_{n}d"
    assert pg.order == 4 * n


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_dihedral_chiral(n):
    pg = dihedral_chiral(n)
    assert pg.name == f"D_{n}"
    assert pg.order == 2 * n


def test_dihedral_invalid_n():
    with pytest.raises(ValueError):
        dihedral(1, "h")


def test_dihedral_chiral_invalid_n():
    with pytest.raises(ValueError):
        dihedral_chiral(1)


def test_dihedral_invalid_variant():
    with pytest.raises(ValueError):
        dihedral(3, "x")


def test_dihedral_h_and_d_share_order():
    """D_nh and D_nd for the same n have the same order but different names."""
    assert dihedral(4, "h").order == dihedral(4, "d").order
    assert dihedral(4, "h") != dihedral(4, "d")


def test_dihedral_chiral_is_half_order_of_full():
    for n in range(2, 7):
        assert dihedral_chiral(n).order == dihedral(n, "h").order // 2


# ---------------------------------------------------------------------------
# Cyclic group factories
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_cyclic_nv(n):
    pg = cyclic(n)
    assert pg.name == f"C_{n}v"
    assert pg.order == 2 * n


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_cyclic_pure(n):
    pg = cyclic_pure(n)
    assert pg.name == f"C_{n}"
    assert pg.order == n


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_cyclic_horizontal(n):
    pg = cyclic_horizontal(n)
    assert pg.name == f"C_{n}h"
    assert pg.order == 2 * n


def test_cyclic_pure_invalid():
    with pytest.raises(ValueError):
        cyclic_pure(0)


def test_cyclic_horizontal_invalid():
    with pytest.raises(ValueError):
        cyclic_horizontal(0)


def test_cyclic_nv_and_nh_share_order():
    """C_nv and C_nh for the same n have equal order but different names."""
    assert cyclic(4).order == cyclic_horizontal(4).order
    assert cyclic(4) != cyclic_horizontal(4)


# ---------------------------------------------------------------------------
# Improper and special groups
# ---------------------------------------------------------------------------


def test_improper_s4():
    pg = improper(2)
    assert pg.name == "S_4"
    assert pg.order == 4


def test_improper_s6():
    pg = improper(3)
    assert pg.name == "S_6"
    assert pg.order == 6


def test_improper_s2_equals_c_i_order():
    """S_2 = C_i: improper(1) and C_I share order 2."""
    pg = improper(1)
    assert pg.name == "S_2"
    assert pg.order == 2
    assert pg.order == C_I.order


def test_improper_invalid():
    with pytest.raises(ValueError):
        improper(0)


def test_c_s_constant():
    assert C_S.name == "C_s"
    assert C_S.order == 2


def test_c_i_constant():
    assert C_I.name == "C_i"
    assert C_I.order == 2


def test_c_s_and_c_i_distinct():
    """C_s and C_i both have order 2 but represent different symmetries."""
    assert C_S != C_I
