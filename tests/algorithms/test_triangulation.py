"""Tests for barycentric subdivision triangulation."""

from __future__ import annotations

import pytest

from polygraph.algorithms.triangulation import (
    BarycentricResult,
    CellOrigin,
    CellType,
    barycentric_subdivision,
    validate_triangulation,
)
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism
from polygraph.generators.johnson import dipyramid, pyramid


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _counts(dm):
    """Return (V, E, F) for a dart map."""
    return (
        len(dm.vertex_orbits()),
        dm.num_edges,
        len(dm.face_orbits()),
    )


# ---------------------------------------------------------------------------
# Basic construction tests
# ---------------------------------------------------------------------------


class TestBarycentricCounts:
    """Verify vertex / edge / face / dart counts after subdivision."""

    @pytest.mark.parametrize(
        "name, make",
        [
            ("tetrahedron", tetrahedron),
            ("cube", cube),
            ("octahedron", octahedron),
            ("dodecahedron", dodecahedron),
            ("icosahedron", icosahedron),
        ],
    )
    def test_platonic_counts(self, name, make):
        dm = make()
        result = barycentric_subdivision(dm)
        sd = result.dart_map

        v, e, f = _counts(dm)
        sv, se, sf = _counts(sd)

        assert sd.num_darts == 6 * dm.num_darts
        assert sv == v + e + f
        assert se == 3 * dm.num_darts  # 6E_orig = 3 * n_orig
        assert sf == 2 * dm.num_darts  # 4E_orig = 2 * n_orig

    @pytest.mark.parametrize("n", [3, 4, 5, 6])
    def test_prism_counts(self, n):
        dm = prism(n)
        result = barycentric_subdivision(dm)
        sd = result.dart_map

        v, e, f = _counts(dm)
        sv, se, sf = _counts(sd)

        assert sv == v + e + f
        assert se == 3 * dm.num_darts
        assert sf == 2 * dm.num_darts

    @pytest.mark.parametrize("n", [3, 4, 5])
    def test_antiprism_counts(self, n):
        dm = antiprism(n)
        result = barycentric_subdivision(dm)
        sd = result.dart_map

        v, e, f = _counts(dm)
        sv, se, sf = _counts(sd)

        assert sv == v + e + f
        assert se == 3 * dm.num_darts
        assert sf == 2 * dm.num_darts

    @pytest.mark.parametrize("n", [3, 4, 5])
    def test_pyramid_counts(self, n):
        dm = pyramid(n)
        result = barycentric_subdivision(dm)
        sd = result.dart_map

        v, e, f = _counts(dm)
        sv, se, sf = _counts(sd)

        assert sv == v + e + f
        assert se == 3 * dm.num_darts
        assert sf == 2 * dm.num_darts


# ---------------------------------------------------------------------------
# All-triangle check
# ---------------------------------------------------------------------------


class TestAllTriangles:
    """Every face in the subdivision must be a triangle."""

    @pytest.mark.parametrize(
        "make",
        [tetrahedron, cube, octahedron, dodecahedron, icosahedron],
    )
    def test_all_faces_triangular(self, make):
        dm = make()
        sd = barycentric_subdivision(dm).dart_map
        for face in sd.face_orbits():
            assert len(face) == 3


# ---------------------------------------------------------------------------
# Topological invariants
# ---------------------------------------------------------------------------


class TestTopology:
    """Euler characteristic and genus must be preserved."""

    @pytest.mark.parametrize(
        "make",
        [
            tetrahedron,
            cube,
            octahedron,
            dodecahedron,
            icosahedron,
            lambda: prism(5),
            lambda: antiprism(4),
            lambda: pyramid(4),
            lambda: dipyramid(5),
        ],
    )
    def test_euler_characteristic_preserved(self, make):
        dm = make()
        sd = barycentric_subdivision(dm).dart_map
        assert sd.euler_characteristic() == dm.euler_characteristic()

    @pytest.mark.parametrize(
        "make",
        [tetrahedron, cube, octahedron, dodecahedron, icosahedron],
    )
    def test_genus_preserved(self, make):
        dm = make()
        sd = barycentric_subdivision(dm).dart_map
        assert sd.genus() == dm.genus() == 0


# ---------------------------------------------------------------------------
# Validation helper
# ---------------------------------------------------------------------------


class TestValidation:
    """The validate_triangulation helper should accept correct subdivisions."""

    @pytest.mark.parametrize(
        "make",
        [tetrahedron, cube, octahedron, dodecahedron, icosahedron],
    )
    def test_validate_passes(self, make):
        dm = make()
        sd = barycentric_subdivision(dm).dart_map
        # Should not raise
        validate_triangulation(dm, sd)


# ---------------------------------------------------------------------------
# Cell map
# ---------------------------------------------------------------------------


class TestCellMap:
    """CellMap must account for all new vertices correctly."""

    @pytest.mark.parametrize(
        "make",
        [tetrahedron, cube, octahedron, dodecahedron, icosahedron],
    )
    def test_cell_type_counts(self, make):
        dm = make()
        result = barycentric_subdivision(dm)

        v, e, f = _counts(dm)
        cm = result.cell_map

        n_vertex = sum(
            1 for c in cm.values() if c.cell_type == CellType.VERTEX
        )
        n_edge = sum(
            1 for c in cm.values() if c.cell_type == CellType.EDGE
        )
        n_face = sum(
            1 for c in cm.values() if c.cell_type == CellType.FACE
        )

        assert n_vertex == v
        assert n_edge == e
        assert n_face == f

    @pytest.mark.parametrize(
        "make",
        [tetrahedron, cube, dodecahedron],
    )
    def test_cell_map_covers_all_vertices(self, make):
        dm = make()
        result = barycentric_subdivision(dm)
        sd = result.dart_map

        # Number of entries in cell_map should equal number of vertices
        assert len(result.cell_map) == len(sd.vertex_orbits())


# ---------------------------------------------------------------------------
# Symmetry preservation
# ---------------------------------------------------------------------------


class TestSymmetryPreservation:
    """Automorphism group order must be preserved by subdivision."""

    @pytest.mark.parametrize(
        "name, make, expected_order",
        [
            ("tetrahedron", tetrahedron, 24),
            ("cube", cube, 48),
            ("octahedron", octahedron, 48),
        ],
    )
    def test_automorphism_group_order_preserved(
        self, name, make, expected_order
    ):
        pynauty = pytest.importorskip("pynauty")
        from polygraph.algorithms.symmetry import (
            automorphism_group_order,
            compute_automorphism_generators,
        )

        dm = make()
        sd = barycentric_subdivision(dm).dart_map

        gens_orig = compute_automorphism_generators(dm)
        order_orig = automorphism_group_order(gens_orig, dm.num_darts)

        gens_sub = compute_automorphism_generators(sd)
        order_sub = automorphism_group_order(gens_sub, sd.num_darts)

        assert order_orig == expected_order
        assert order_sub == order_orig


# ---------------------------------------------------------------------------
# Specific known values
# ---------------------------------------------------------------------------


class TestKnownValues:
    """Spot-check concrete numbers for well-known polyhedra."""

    def test_tetrahedron_subdivision(self):
        dm = tetrahedron()
        sd = barycentric_subdivision(dm).dart_map
        # Tetrahedron: V=4, E=6, F=4, n=12
        assert _counts(sd) == (4 + 6 + 4, 3 * 12, 2 * 12)  # (14, 36, 24)

    def test_cube_subdivision(self):
        dm = cube()
        sd = barycentric_subdivision(dm).dart_map
        # Cube: V=8, E=12, F=6, n=24
        assert _counts(sd) == (8 + 12 + 6, 3 * 24, 2 * 24)  # (26, 72, 48)

    def test_dodecahedron_subdivision(self):
        dm = dodecahedron()
        sd = barycentric_subdivision(dm).dart_map
        # Dodecahedron: V=20, E=30, F=12, n=60
        assert _counts(sd) == (20 + 30 + 12, 3 * 60, 2 * 60)  # (62, 180, 120)
