"""Tests for DartMap automorphism computation and orbit analysis.

Validation targets from the roadmap:
  Tetrahedron  : |Aut| = 24,  1 vertex orbit, 1 edge orbit, 1 face orbit
  Cube         : |Aut| = 48,  1 vertex orbit, 1 edge orbit, 1 face orbit
  Icosahedron  : |Aut| = 120, 1 vertex orbit, 1 edge orbit, 1 face orbit
  Prism(5)     : |Aut| = 20,  2 vertex orbits, 2 edge orbits, 2 face orbits
"""

from __future__ import annotations

import pytest

from polygraph.algorithms.symmetry import (
    automorphism_group_order,
    compute_automorphism_generators,
    dart_orbits,
    edge_orbits,
    face_orbits,
    is_orientation_preserving,
    vertex_orbits,
)
from polygraph.generators.platonic import (
    cube,
    icosahedron,
    tetrahedron,
)
from polygraph.generators.prisms import prism
from polygraph.structures.permutation import Permutation

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tet():
    return tetrahedron()


@pytest.fixture
def cub():
    return cube()


@pytest.fixture
def ico():
    return icosahedron()


@pytest.fixture
def prism5():
    return prism(5)


# ---------------------------------------------------------------------------
# Generator validity
# ---------------------------------------------------------------------------

def test_generators_are_valid_permutations(tet):
    gens = compute_automorphism_generators(tet)
    assert all(isinstance(g, Permutation) for g in gens)
    assert all(len(g) == tet.num_darts for g in gens)


def test_generators_commute_with_phi_or_reverse_phi(tet):
    """Full-group generators either preserve or reverse face traversal."""
    gens = compute_automorphism_generators(tet)
    phi = [tet._sigma_inv[tet.alpha[d]] for d in range(tet.num_darts)]
    for g in gens:
        preserves = all(g[phi[d]] == phi[g[d]] for d in range(tet.num_darts))
        reverses = all(
            phi[g[phi[d]]] == g[d] for d in range(tet.num_darts)
        )
        assert preserves or reverses, (
            "Generator neither commutes with phi nor reverses phi"
        )


def test_generators_commute_with_alpha(tet):
    gens = compute_automorphism_generators(tet)
    for g in gens:
        for d in range(tet.num_darts):
            assert g[tet.alpha[d]] == tet.alpha[g[d]], (
                f"Generator does not commute with alpha at dart {d}"
            )


# ---------------------------------------------------------------------------
# Group order: tetrahedron
# ---------------------------------------------------------------------------

def test_tetrahedron_group_order_from_pynauty(tet):
    _, order = __import__(
        "polygraph.interop.bliss_adapter",
        fromlist=["dartmap_automorphisms"],
    ).dartmap_automorphisms(tet)
    assert order == 24


def test_tetrahedron_group_order_from_generators(tet):
    gens = compute_automorphism_generators(tet)
    assert automorphism_group_order(gens, tet.num_darts) == 24


# ---------------------------------------------------------------------------
# Group order: cube
# ---------------------------------------------------------------------------

def test_cube_group_order(cub):
    _, order = __import__(
        "polygraph.interop.bliss_adapter",
        fromlist=["dartmap_automorphisms"],
    ).dartmap_automorphisms(cub)
    assert order == 48


# ---------------------------------------------------------------------------
# Group order: icosahedron
# ---------------------------------------------------------------------------

def test_icosahedron_group_order(ico):
    _, order = __import__(
        "polygraph.interop.bliss_adapter",
        fromlist=["dartmap_automorphisms"],
    ).dartmap_automorphisms(ico)
    assert order == 120


# ---------------------------------------------------------------------------
# Orbit counts: Platonic solids have 1 orbit each
# ---------------------------------------------------------------------------

def test_tetrahedron_single_vertex_orbit(tet):
    gens = compute_automorphism_generators(tet)
    assert len(vertex_orbits(gens, tet)) == 1


def test_tetrahedron_single_edge_orbit(tet):
    gens = compute_automorphism_generators(tet)
    assert len(edge_orbits(gens, tet)) == 1


def test_tetrahedron_single_face_orbit(tet):
    gens = compute_automorphism_generators(tet)
    assert len(face_orbits(gens, tet)) == 1


def test_cube_single_vertex_orbit(cub):
    gens = compute_automorphism_generators(cub)
    assert len(vertex_orbits(gens, cub)) == 1


def test_cube_single_edge_orbit(cub):
    gens = compute_automorphism_generators(cub)
    assert len(edge_orbits(gens, cub)) == 1


def test_cube_single_face_orbit(cub):
    gens = compute_automorphism_generators(cub)
    assert len(face_orbits(gens, cub)) == 1


def test_icosahedron_single_vertex_orbit(ico):
    gens = compute_automorphism_generators(ico)
    assert len(vertex_orbits(gens, ico)) == 1


def test_icosahedron_single_edge_orbit(ico):
    gens = compute_automorphism_generators(ico)
    assert len(edge_orbits(gens, ico)) == 1


def test_icosahedron_single_face_orbit(ico):
    gens = compute_automorphism_generators(ico)
    assert len(face_orbits(gens, ico)) == 1


# ---------------------------------------------------------------------------
# Orbit counts: prism(5) has 2 orbits each
# ---------------------------------------------------------------------------

def test_prism5_group_order(prism5):
    _, order = __import__(
        "polygraph.interop.bliss_adapter",
        fromlist=["dartmap_automorphisms"],
    ).dartmap_automorphisms(prism5)
    assert order == 20


def test_prism5_one_vertex_orbit(prism5):
    # Full D_5h symmetry group (order 20) contains the horizontal C_2
    # rotations that map top vertices to bottom vertices, giving a single
    # vertex orbit.  (The roadmap's "2 vertex orbits" is incorrect for the
    # full group; it only holds for the orientation-preserving D_5 subgroup
    # which does NOT include sigma-reversing reflections.)
    gens = compute_automorphism_generators(prism5)
    assert len(vertex_orbits(gens, prism5)) == 1


def test_prism5_two_edge_orbits(prism5):
    gens = compute_automorphism_generators(prism5)
    assert len(edge_orbits(gens, prism5)) == 2


def test_prism5_two_face_orbits(prism5):
    gens = compute_automorphism_generators(prism5)
    assert len(face_orbits(gens, prism5)) == 2


# ---------------------------------------------------------------------------
# Orbit partition sanity checks
# ---------------------------------------------------------------------------

def test_vertex_orbits_cover_all_vertices(tet):
    gens = compute_automorphism_generators(tet)
    orbits = vertex_orbits(gens, tet)
    # Each orbit element is a vertex representative.
    all_reps = [d for orbit in orbits for d in orbit]
    # Count: should equal number of vertices in the tetrahedron (4).
    assert len(all_reps) == 4
    assert len(set(all_reps)) == 4  # no duplicates


def test_edge_orbits_cover_all_edges(tet):
    gens = compute_automorphism_generators(tet)
    orbits = edge_orbits(gens, tet)
    all_reps = [d for orbit in orbits for d in orbit]
    assert len(all_reps) == tet.num_edges
    assert len(set(all_reps)) == tet.num_edges


def test_face_orbits_cover_all_faces(tet):
    gens = compute_automorphism_generators(tet)
    orbits = face_orbits(gens, tet)
    all_reps = [d for orbit in orbits for d in orbit]
    num_faces = len(tet.face_orbits())
    assert len(all_reps) == num_faces
    assert len(set(all_reps)) == num_faces


def test_dart_orbits_cover_all_darts(tet):
    gens = compute_automorphism_generators(tet)
    orbits = dart_orbits(gens, tet)
    all_darts = sorted(d for orbit in orbits for d in orbit)
    assert all_darts == list(range(tet.num_darts))


# ---------------------------------------------------------------------------
# is_orientation_preserving
# ---------------------------------------------------------------------------

def test_some_generators_are_orientation_reversing(tet):
    """The full tetrahedral group T_d contains reflections (order 24 > 12)."""
    gens = compute_automorphism_generators(tet)
    reversing = [g for g in gens if not is_orientation_preserving(g, tet)]
    assert len(reversing) > 0, (
        "Expected at least one orientation-reversing generator for T_d"
    )


def test_identity_is_orientation_preserving(tet):
    identity = Permutation.identity(tet.num_darts)
    assert is_orientation_preserving(identity, tet)


def test_is_orientation_preserving_rejects_wrong_size(tet):
    wrong = Permutation.identity(tet.num_darts + 2)
    with pytest.raises(ValueError, match="domain size"):
        is_orientation_preserving(wrong, tet)


# ---------------------------------------------------------------------------
# automorphism_group_order edge cases
# ---------------------------------------------------------------------------

def test_group_order_empty_generators_is_one(tet):
    assert automorphism_group_order([], tet.num_darts) == 1
