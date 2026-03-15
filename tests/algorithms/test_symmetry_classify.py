"""Tests for symmetry point group classification."""

import pytest

pytest.importorskip("pynauty")

from polygraph.algorithms.symmetry import (  # noqa: E402
    classify_symmetry,
    compute_automorphism_generators,
)
from polygraph.algorithms.symmetry.classify import UnknownSymmetryError
from polygraph.algorithms.symmetry.point_groups import PointGroup
from polygraph.generators.johnson import dipyramid, pyramid
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism

# ---------------------------------------------------------------------------
# Platonic solids
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_name, expected_order",
    [
        (tetrahedron, "T_d", 24),
        (cube, "O_h", 48),
        (octahedron, "O_h", 48),
        (dodecahedron, "I_h", 120),
        (icosahedron, "I_h", 120),
    ],
)
def test_platonic_classify(factory, expected_name, expected_order):
    dm = factory()
    gens = compute_automorphism_generators(dm)
    pg = classify_symmetry(gens, dm)
    assert pg.name == expected_name
    assert pg.order == expected_order


# ---------------------------------------------------------------------------
# Prisms: D_nh
# n=4 excluded — prism(4) == cube, classified as O_h
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [3, 5, 6, 7])
def test_prism_classify(n):
    dm = prism(n)
    gens = compute_automorphism_generators(dm)
    pg = classify_symmetry(gens, dm)
    assert pg.name == f"D_{n}h"
    assert pg.order == 4 * n


# ---------------------------------------------------------------------------
# Antiprisms: D_nd
# n=3 excluded — antiprism(3) == octahedron, classified as O_h
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [4, 5, 6])
def test_antiprism_classify(n):
    dm = antiprism(n)
    gens = compute_automorphism_generators(dm)
    pg = classify_symmetry(gens, dm)
    assert pg.name == f"D_{n}d"
    assert pg.order == 4 * n


# ---------------------------------------------------------------------------
# Dipyramids: D_nh
# n=4 excluded — dipyramid(4) == octahedron, classified as O_h
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [3, 5, 6])
def test_dipyramid_classify(n):
    dm = dipyramid(n)
    gens = compute_automorphism_generators(dm)
    pg = classify_symmetry(gens, dm)
    assert pg.name == f"D_{n}h"
    assert pg.order == 4 * n


# ---------------------------------------------------------------------------
# Pyramids: C_nv
# n=3 excluded — pyramid(3) == tetrahedron, classified as T_d
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [4, 5, 6])
def test_pyramid_classify(n):
    dm = pyramid(n)
    gens = compute_automorphism_generators(dm)
    pg = classify_symmetry(gens, dm)
    assert pg.name == f"C_{n}v"
    assert pg.order == 2 * n


# ---------------------------------------------------------------------------
# Structural / API tests
# ---------------------------------------------------------------------------


def test_classify_returns_point_group_instance():
    dm = tetrahedron()
    gens = compute_automorphism_generators(dm)
    pg = classify_symmetry(gens, dm)
    assert isinstance(pg, PointGroup)


def test_classify_accepts_precomputed_order():
    """Passing order= bypasses BFS closure; result must match."""
    dm = cube()
    gens = compute_automorphism_generators(dm)
    pg_auto = classify_symmetry(gens, dm)
    pg_manual = classify_symmetry(gens, dm, order=48)
    assert pg_auto == pg_manual
    assert pg_manual.name == "O_h"


# ---------------------------------------------------------------------------
# Regression tests for lateral-face disambiguation
# ---------------------------------------------------------------------------


def test_prism3_is_d3h_not_d3d():
    """prism(3) has triangular caps and quad lateral faces → D_3h."""
    dm = prism(3)
    gens = compute_automorphism_generators(dm)
    assert classify_symmetry(gens, dm).name == "D_3h"


def test_antiprism4_is_d4d_not_d4h():
    """antiprism(4) has square caps and triangular lateral faces → D_4d."""
    dm = antiprism(4)
    gens = compute_automorphism_generators(dm)
    assert classify_symmetry(gens, dm).name == "D_4d"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


def test_empty_generators_raises_unknown_symmetry_error():
    """Trivial automorphism group (order 1) has no matching point group."""
    dm = tetrahedron()
    with pytest.raises(UnknownSymmetryError):
        classify_symmetry([], dm)
