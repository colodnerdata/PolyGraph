"""Tests for topology traversal helpers."""

from __future__ import annotations

import pytest

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation
from polygraph.structures.traversal import (
    adjacent_face_pairs,
    all_edge_orbits,
    all_face_orbits,
    all_vertex_orbits,
    edge_darts,
    face_darts,
    faces_incident_to_vertex,
    neighbors,
    orbit,
    vertex_darts,
    vertices_of_face,
)


def _tetrahedron_map() -> DartMap:
    faces = [
        [0, 1, 2],
        [0, 3, 1],
        [1, 3, 2],
        [0, 2, 3],
    ]
    return DartMap.from_face_lists(faces, num_vertices=4)


def test_orbit_yields_full_cycle_without_repeating_start() -> None:
    perm = Permutation.from_sequence([1, 2, 0, 4, 3])
    assert list(orbit(0, perm)) == [0, 1, 2]
    assert list(orbit(3, perm)) == [3, 4]


def test_orbit_raises_for_invalid_start_index() -> None:
    perm = Permutation.identity(3)
    with pytest.raises(IndexError):
        list(orbit(-1, perm))
    with pytest.raises(IndexError):
        list(orbit(3, perm))


def test_permutation_orbit_raises_for_invalid_start() -> None:
    perm = Permutation.from_sequence([1, 2, 0])
    with pytest.raises(IndexError):
        list(perm.orbit(-1))
    with pytest.raises(IndexError):
        list(perm.orbit(3))


def test_permutation_identity_raises_for_negative_n() -> None:
    with pytest.raises(ValueError):
        Permutation.identity(-1)


def test_traversal_functions_raise_for_out_of_range_dart() -> None:
    dm = _tetrahedron_map()
    with pytest.raises(IndexError):
        list(vertex_darts(dm, -1))
    with pytest.raises(IndexError):
        list(vertex_darts(dm, dm.num_darts))
    with pytest.raises(IndexError):
        list(face_darts(dm, -1))
    with pytest.raises(IndexError):
        list(face_darts(dm, dm.num_darts))
    with pytest.raises(IndexError):
        list(edge_darts(dm, -1))
    with pytest.raises(IndexError):
        list(edge_darts(dm, dm.num_darts))


def test_vertex_face_and_edge_darts_are_consistent() -> None:
    dm = _tetrahedron_map()
    v_rep = next(all_vertex_orbits(dm))
    f_rep = next(all_face_orbits(dm))
    e_rep = next(all_edge_orbits(dm))

    v_cycle = list(vertex_darts(dm, v_rep))
    f_cycle = list(face_darts(dm, f_rep))
    e_cycle = list(edge_darts(dm, e_rep))

    assert len(v_cycle) == 3
    assert dm.sigma[v_cycle[-1]] == v_cycle[0]

    assert len(f_cycle) == 3
    assert dm.phi(f_cycle[-1]) == f_cycle[0]

    assert e_cycle == [e_rep, dm.alpha[e_rep]]
    assert dm.alpha[e_cycle[1]] == e_cycle[0]


def test_orbit_enumerators_match_dart_map_counts_and_cover_darts() -> None:
    dm = _tetrahedron_map()
    v_reps = list(all_vertex_orbits(dm))
    f_reps = list(all_face_orbits(dm))
    e_reps = list(all_edge_orbits(dm))

    assert len(v_reps) == len(dm.vertex_orbits())
    assert len(f_reps) == len(dm.face_orbits())
    assert len(e_reps) == dm.num_edges

    vertex_covered = sorted(d for r in v_reps for d in vertex_darts(dm, r))
    face_covered = sorted(d for r in f_reps for d in face_darts(dm, r))
    assert vertex_covered == list(range(dm.num_darts))
    assert face_covered == list(range(dm.num_darts))

    edge_pairs = [tuple(sorted(edge_darts(dm, r))) for r in e_reps]
    assert len(edge_pairs) == len(set(edge_pairs)) == dm.num_edges
    assert all(r < dm.alpha[r] for r in e_reps)


def test_neighbors_preserves_cyclic_order_and_multiplicity() -> None:
    dm = _tetrahedron_map()
    v_rep = next(all_vertex_orbits(dm))
    expected = [dm.alpha[d] for d in vertex_darts(dm, v_rep)]
    observed = list(neighbors(dm, v_rep))

    assert observed == expected
    assert len(observed) == len(list(vertex_darts(dm, v_rep)))


def test_skeleton_traversal_methods_align_with_rep_dart_model() -> None:
    dm = _tetrahedron_map()
    v_rep = next(all_vertex_orbits(dm))
    f_rep = next(all_face_orbits(dm))

    assert list(vertices_of_face(dm, f_rep)) == list(face_darts(dm, f_rep))
    assert list(faces_incident_to_vertex(dm, v_rep)) == list(
        vertex_darts(dm, v_rep)
    )

    pairs = list(adjacent_face_pairs(dm))
    assert len(pairs) == dm.num_edges
    assert [a for a, _ in pairs] == list(all_edge_orbits(dm))
    for a, b in pairs:
        assert dm.alpha[a] == b
        assert dm.alpha[b] == a
