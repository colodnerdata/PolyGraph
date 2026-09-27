"""Tests for the symmetry automorphism layer."""

import pytest

pytest.importorskip("pynauty")

from polygraph.algorithms.symmetry import (
    automorphism_group_order,
    compute_automorphism_generators,
    dart_orbits,
    edge_orbits,
    face_orbits,
    is_orientation_preserving,
    vertex_orbits,
)
from polygraph.generators.johnson import (
    dipyramid,
    gyroelongated_square_bipyramid,
    pyramid,
    snub_disphenoid,
    triaugmented_triangular_prism,
)
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import prism
from polygraph.structures.permutation import Permutation

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _phi_perm(dm):
    """Return the full phi permutation as a list."""
    return [dm.phi(d) for d in range(dm.num_darts)]


def _phi_inv_perm(dm):
    """Return the inverse phi permutation as a list."""
    phi = _phi_perm(dm)
    inv = [0] * dm.num_darts
    for d, p in enumerate(phi):
        inv[p] = d
    return inv


# ---------------------------------------------------------------------------
# Generator validity
# ---------------------------------------------------------------------------


def _assert_generator_valid(gen, dm):
    """Assert *gen* is a valid dart-map automorphism."""
    alpha = dm.alpha
    phi = _phi_perm(dm)
    phi_inv = _phi_inv_perm(dm)

    # Must commute with alpha.
    for d in range(dm.num_darts):
        assert gen[alpha[d]] == alpha[gen[d]], (
            f"Generator does not commute with alpha at dart {d}"
        )

    # Must commute with phi OR reverse phi (one of the two must hold for
    # every dart simultaneously).
    preserving = all(gen[phi[d]] == phi[gen[d]] for d in range(dm.num_darts))
    reversing = all(
        gen[phi[d]] == phi_inv[gen[d]] for d in range(dm.num_darts)
    )
    assert preserving or reversing, (
        "Generator neither commutes with phi nor reverses phi"
    )


def test_tetrahedron_generators_valid():
    dm = tetrahedron()
    for gen in compute_automorphism_generators(dm):
        _assert_generator_valid(gen, dm)


def test_cube_generators_valid():
    dm = cube()
    for gen in compute_automorphism_generators(dm):
        _assert_generator_valid(gen, dm)


def test_icosahedron_generators_valid():
    dm = icosahedron()
    for gen in compute_automorphism_generators(dm):
        _assert_generator_valid(gen, dm)


# ---------------------------------------------------------------------------
# Group orders
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_order",
    [
        (tetrahedron, 24),
        (cube, 48),
        (octahedron, 48),
        (dodecahedron, 120),
        (icosahedron, 120),
    ],
)
def test_platonic_group_order(factory, expected_order):
    dm = factory()
    _gens, order = __import__(
        "polygraph.interop.pynauty_adapter", fromlist=["dartmap_automorphisms"]
    ).dartmap_automorphisms(dm)
    assert order == expected_order


def test_bfs_group_order_matches_pynauty():
    from polygraph.interop.pynauty_adapter import dartmap_automorphisms

    dm = tetrahedron()
    gens, pynauty_order = dartmap_automorphisms(dm)
    bfs_order = automorphism_group_order(gens, dm.num_darts)
    assert bfs_order == pynauty_order == 24


def test_empty_generators_order_one():
    dm = tetrahedron()
    assert automorphism_group_order([], dm.num_darts) == 1


# ---------------------------------------------------------------------------
# Orbit counts — Platonic solids (all 1 each)
# ---------------------------------------------------------------------------


_PLATONICS = [tetrahedron, cube, octahedron, dodecahedron, icosahedron]


@pytest.mark.parametrize("factory", _PLATONICS)
def test_platonic_one_vertex_orbit(factory):
    dm = factory()
    gens = compute_automorphism_generators(dm)
    assert len(vertex_orbits(gens, dm)) == 1


@pytest.mark.parametrize("factory", _PLATONICS)
def test_platonic_one_edge_orbit(factory):
    dm = factory()
    gens = compute_automorphism_generators(dm)
    assert len(edge_orbits(gens, dm)) == 1


@pytest.mark.parametrize("factory", _PLATONICS)
def test_platonic_one_face_orbit(factory):
    dm = factory()
    gens = compute_automorphism_generators(dm)
    assert len(face_orbits(gens, dm)) == 1


# ---------------------------------------------------------------------------
# Orbit counts — triangular prism
# ---------------------------------------------------------------------------


def test_prism3_vertex_orbits():
    dm = prism(3)
    gens = compute_automorphism_generators(dm)
    assert len(vertex_orbits(gens, dm)) == 1


def test_prism3_edge_orbits():
    dm = prism(3)
    gens = compute_automorphism_generators(dm)
    assert len(edge_orbits(gens, dm)) == 2


def test_prism3_face_orbits():
    dm = prism(3)
    gens = compute_automorphism_generators(dm)
    assert len(face_orbits(gens, dm)) == 2


# ---------------------------------------------------------------------------
# Orbit counts — non-vertex-transitive solids
#
# Regression tests for orientation-reversing generators.  Under the pynauty
# phi-convention a reflection maps sigma-orbits to alpha-images of
# sigma-orbits, so a naive vertex action merges distinct vertex orbits.
# Every solid below has a reflection and vertices of two distinct degrees,
# so a wrong action collapses the count to 1.
# ---------------------------------------------------------------------------


def _degree_of_rep(dm, rep):
    """Return the degree of the vertex whose representative dart is *rep*."""
    for cycle in dm.vertex_orbits():
        if min(cycle) == rep:
            return len(cycle)
    raise AssertionError(f"{rep} is not a vertex representative")


@pytest.mark.parametrize(
    "factory, expected_vertex_orbits",
    [
        (lambda: pyramid(4), 2),
        (lambda: pyramid(5), 2),
        (lambda: dipyramid(3), 2),
        (lambda: dipyramid(5), 2),
        (snub_disphenoid, 2),
        (triaugmented_triangular_prism, 2),
        (gyroelongated_square_bipyramid, 2),
    ],
    ids=["pyramid4", "pyramid5", "dipyramid3", "dipyramid5", "J84", "J51",
         "J17"],
)
def test_vertex_orbit_count_with_reflections(
    factory, expected_vertex_orbits
):
    dm = factory()
    gens = compute_automorphism_generators(dm)
    assert any(not is_orientation_preserving(g, dm) for g in gens)
    assert len(vertex_orbits(gens, dm)) == expected_vertex_orbits


@pytest.mark.parametrize(
    "factory",
    [lambda: pyramid(4), lambda: dipyramid(3), snub_disphenoid],
    ids=["pyramid4", "dipyramid3", "J84"],
)
def test_vertex_orbits_have_uniform_degree(factory):
    dm = factory()
    gens = compute_automorphism_generators(dm)
    for orbit in vertex_orbits(gens, dm):
        degrees = {_degree_of_rep(dm, rep) for rep in orbit}
        assert len(degrees) == 1, f"orbit {orbit} mixes degrees {degrees}"


def test_reversing_generator_twisted_by_alpha_preserves_vertices():
    """Check that ``alpha ∘ g`` maps sigma-orbits onto sigma-orbits."""
    dm = dipyramid(3)
    gens = compute_automorphism_generators(dm)
    reversing = [g for g in gens if not is_orientation_preserving(g, dm)]
    assert reversing

    dart_to_vertex = [0] * dm.num_darts
    for cycle in dm.vertex_orbits():
        for d in cycle:
            dart_to_vertex[d] = min(cycle)

    for g in reversing:
        for cycle in dm.vertex_orbits():
            images = {dart_to_vertex[dm.alpha[g[d]]] for d in cycle}
            assert len(images) == 1


# ---------------------------------------------------------------------------
# Orbit coverage
# ---------------------------------------------------------------------------


def test_vertex_orbits_cover_all_vertices():
    dm = cube()
    gens = compute_automorphism_generators(dm)
    orbits = vertex_orbits(gens, dm)
    all_reps = {min(cycle) for cycle in dm.vertex_orbits()}
    covered = {rep for orbit in orbits for rep in orbit}
    assert covered == all_reps


def test_edge_orbits_cover_all_edges():
    dm = cube()
    gens = compute_automorphism_generators(dm)
    orbits = edge_orbits(gens, dm)
    all_edges = {d for d in range(dm.num_darts) if d < dm.alpha[d]}
    covered = {rep for orbit in orbits for rep in orbit}
    assert covered == all_edges


def test_face_orbits_cover_all_faces():
    dm = cube()
    gens = compute_automorphism_generators(dm)
    orbits = face_orbits(gens, dm)
    all_reps = {min(cycle) for cycle in dm.face_orbits()}
    covered = {rep for orbit in orbits for rep in orbit}
    assert covered == all_reps


# ---------------------------------------------------------------------------
# Dart orbits
# ---------------------------------------------------------------------------


def test_dart_orbits_cover_all_darts():
    dm = tetrahedron()
    gens = compute_automorphism_generators(dm)
    orbits = dart_orbits(gens, dm)
    covered = {d for orbit in orbits for d in orbit}
    assert covered == set(range(dm.num_darts))


@pytest.mark.parametrize("factory", _PLATONICS)
def test_platonic_one_dart_orbit(factory):
    # All Platonic solids are flag-transitive: all darts are in one orbit.
    dm = factory()
    gens = compute_automorphism_generators(dm)
    assert len(dart_orbits(gens, dm)) == 1


# ---------------------------------------------------------------------------
# Orientation
# ---------------------------------------------------------------------------


def test_identity_is_orientation_preserving():
    dm = tetrahedron()
    identity = Permutation.from_sequence(range(dm.num_darts))
    assert is_orientation_preserving(identity, dm)


def test_some_generators_are_orientation_reversing():
    # Mirror-symmetric polyhedra must have at least one reflection generator.
    dm = tetrahedron()
    gens = compute_automorphism_generators(dm)
    assert any(not is_orientation_preserving(g, dm) for g in gens), (
        "All generators are orientation-preserving — reflections are missing"
    )
