"""Tests for symmetry detection and orbit computations."""

from __future__ import annotations

import pytest

from polygraph.algorithms.symmetry.automorphisms import (
    automorphism_group_order,
    compute_automorphism_generators,
)
from polygraph.algorithms.symmetry.orbits import (
    edge_orbits,
    face_orbits,
    vertex_orbits,
)
from polygraph.generators.platonic import cube, icosahedron, tetrahedron
from polygraph.generators.prisms import prism


@pytest.mark.parametrize(
    (
        "factory",
        "expected_order",
        "expected_vertex_orbits",
        "expected_edge_orbits",
        "expected_face_orbits",
    ),
    [
        (tetrahedron, 12, 1, 1, 1),
        (cube, 24, 1, 1, 1),
        (icosahedron, 60, 1, 1, 1),
        (lambda: prism(5), 10, 1, 2, 2),
    ],
)
def test_symmetry_detection_matches_known_group_orders(
    factory,
    expected_order: int,
    expected_vertex_orbits: int,
    expected_edge_orbits: int,
    expected_face_orbits: int,
) -> None:
    dm = factory()
    generators = compute_automorphism_generators(dm)

    assert automorphism_group_order(generators, dm.num_darts) == expected_order
    assert len(vertex_orbits(generators, dm)) == expected_vertex_orbits
    assert len(edge_orbits(generators, dm)) == expected_edge_orbits
    assert len(face_orbits(generators, dm)) == expected_face_orbits
