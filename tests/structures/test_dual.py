"""Tests for dual-map construction and dual isomorphism behavior."""

from __future__ import annotations

import pytest

from polygraph.generators.johnson import dipyramid, pyramid
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism
from polygraph.structures.dual import dual_of, is_isomorphism


@pytest.mark.parametrize(
    ("factory", "expected"),
    [
        (tetrahedron, (4, 6, 4)),
        (cube, (6, 12, 8)),
        (octahedron, (8, 12, 6)),
        (dodecahedron, (12, 30, 20)),
        (icosahedron, (20, 30, 12)),
    ],
)
def test_dual_counts_match_expected(factory, expected) -> None:
    dm = factory()
    ddual = dual_of(dm)

    assert len(ddual.vertex_orbits()) == expected[0]
    assert ddual.num_edges == expected[1]
    assert len(ddual.face_orbits()) == expected[2]


def test_dual_sigma_is_phi_and_alpha_unchanged(platonic_dm) -> None:
    ddual = dual_of(platonic_dm)

    assert ddual.alpha == platonic_dm.alpha
    assert [ddual.sigma[d] for d in range(platonic_dm.num_darts)] == [
        platonic_dm.phi(d) for d in range(platonic_dm.num_darts)
    ]


def test_dual_preserves_topological_invariants(platonic_dm) -> None:
    ddual = dual_of(platonic_dm)

    assert ddual.euler_characteristic() == 2
    assert ddual.genus() == 0
    assert ddual.num_darts == platonic_dm.num_darts


def test_dual_vertex_degrees_match_original_face_sizes(platonic_dm) -> None:
    ddual = dual_of(platonic_dm)

    dual_vertex_degrees = sorted(len(orbit) for orbit in ddual.vertex_orbits())
    orig_face_sizes = sorted(len(orbit) for orbit in platonic_dm.face_orbits())

    dual_face_sizes = sorted(len(orbit) for orbit in ddual.face_orbits())
    orig_vertex_degrees = sorted(
        len(orbit) for orbit in platonic_dm.vertex_orbits()
    )

    assert dual_vertex_degrees == orig_face_sizes
    assert dual_face_sizes == orig_vertex_degrees


@pytest.mark.parametrize(
    "dm",
    [
        tetrahedron(),
        cube(),
        octahedron(),
        dodecahedron(),
        icosahedron(),
        prism(3),
        prism(5),
        antiprism(3),
        antiprism(5),
        pyramid(3),
        pyramid(6),
        dipyramid(3),
        dipyramid(6),
    ],
)
def test_double_dual_is_isomorphic_via_alpha(dm) -> None:
    ddm = dual_of(dual_of(dm))

    assert len(ddm.vertex_orbits()) == len(dm.vertex_orbits())
    assert ddm.num_edges == dm.num_edges
    assert len(ddm.face_orbits()) == len(dm.face_orbits())
    assert is_isomorphism(dm, ddm, dm.alpha)


def test_identity_is_isomorphism_for_any_map(platonic_dm) -> None:
    identity = list(range(platonic_dm.num_darts))
    assert is_isomorphism(platonic_dm, platonic_dm, identity)


def test_is_isomorphism_rejects_wrong_length() -> None:
    dm = cube()
    assert not is_isomorphism(dm, dm, [0, 1, 2])


def test_is_isomorphism_rejects_non_permutation() -> None:
    dm = tetrahedron()
    not_a_perm = [0] * dm.num_darts
    assert not is_isomorphism(dm, dm, not_a_perm)
