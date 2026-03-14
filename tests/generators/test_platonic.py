"""Tests for Platonic solid generators."""

from __future__ import annotations

import pytest

from polygraph.generators.platonic import dodecahedron, icosahedron, tetrahedron
from polygraph.structures.dart_map import DartMap


@pytest.mark.parametrize(
    ("generator", "expected_counts", "expected_face_size", "expected_vertex_degree"),
    [
        (tetrahedron, (4, 6, 4), 3, 3),
        (dodecahedron, (20, 30, 12), 5, 3),
        (icosahedron, (12, 30, 20), 3, 5),
    ],
)
def test_platonic_generators_return_expected_topology(
    generator,
    expected_counts: tuple[int, int, int],
    expected_face_size: int,
    expected_vertex_degree: int,
) -> None:
    dm = generator()

    assert isinstance(dm, DartMap)

    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())
    assert (num_vertices, num_edges, num_faces) == expected_counts

    face_sizes = {len(face) for face in dm.face_orbits()}
    assert face_sizes == {expected_face_size}

    vertex_degrees = {len(vertex) for vertex in dm.vertex_orbits()}
    assert vertex_degrees == {expected_vertex_degree}

    assert dm.euler_characteristic() == 2


def test_platonic_generators_pair_edges_in_opposite_directions() -> None:
    for generator in (tetrahedron, dodecahedron, icosahedron):
        dm = generator()

        # alpha involution guarantees opposite orientation pairing for each edge.
        for dart in range(dm.num_darts):
            twin = dm.alpha[dart]
            assert twin != dart
            assert dm.alpha[twin] == dart

        # each canonical undirected edge should be represented exactly once.
        edge_ids = {dm.edge_of_dart(dart) for dart in range(dm.num_darts)}
        assert len(edge_ids) == dm.num_edges
