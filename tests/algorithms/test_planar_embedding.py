"""Tests for Phase 7: planar embedding infrastructure."""

from __future__ import annotations

import pytest

from polygraph.algorithms.planar.embedding import PlanarEmbeddingView
from polygraph.algorithms.planar.outer_face import (
    choose_outer_face,
    outer_face_anchors,
)
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism
from polygraph.structures.dart_map import DartMap

# ---------------------------------------------------------------------------
# Shared fixture: square pyramid with known geometry
#   base: 0=(0,0), 1=(2,0), 2=(2,2), 3=(0,2); apex 4=(1,1)
#   interior triangles wound CCW; outer face [0,3,2,1] wound CW
# ---------------------------------------------------------------------------

PYRAMID_FACES = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4], [0, 3, 2, 1]]


@pytest.fixture(scope="module")
def pyramid_dm():
    return DartMap.from_face_lists(PYRAMID_FACES, 5)


@pytest.fixture(scope="module")
def pyramid_view(pyramid_dm):
    return PlanarEmbeddingView(pyramid_dm)


# ---------------------------------------------------------------------------
# PlanarEmbeddingView — counts
# ---------------------------------------------------------------------------


class TestPlanarEmbeddingViewCounts:
    @pytest.mark.parametrize(
        "make, v, e, f",
        [
            (tetrahedron, 4, 6, 4),
            (cube, 8, 12, 6),
            (octahedron, 6, 12, 8),
            (dodecahedron, 20, 30, 12),
            (icosahedron, 12, 30, 20),
        ],
    )
    def test_platonic_counts(self, make, v, e, f):
        view = PlanarEmbeddingView(make())
        assert view.num_vertices == v
        assert view.num_edges == e
        assert view.num_faces == f

    def test_pyramid_counts(self, pyramid_view):
        assert pyramid_view.num_vertices == 5
        assert pyramid_view.num_edges == 8
        assert pyramid_view.num_faces == 5


# ---------------------------------------------------------------------------
# PlanarEmbeddingView — ordered_neighbors
# ---------------------------------------------------------------------------


class TestOrderedNeighbors:
    def test_sigma_orbit_order_matches(self, pyramid_view, pyramid_dm):
        """ordered_neighbors must agree with the sigma orbit traversal."""
        for v in range(pyramid_view.num_vertices):
            rep = pyramid_view._vertex_to_rep[v]
            sigma_neighbors = []
            d = rep
            while True:
                sigma_neighbors.append(
                    pyramid_view._dart_to_vertex[pyramid_dm.alpha[d]]
                )
                d = pyramid_dm.sigma[d]
                if d == rep:
                    break
            assert pyramid_view.ordered_neighbors(v) == sigma_neighbors

    def test_apex_ccw_order(self, pyramid_view):
        """Apex neighbours should appear in CCW angular order.

        In PYRAMID_FACES, vertex label 4 is the apex.  After vertex_orbits()
        scan, the apex receives vertex ID 2 (first dart of its sigma orbit is
        dart 2).  From apex (1,1): labels 0→225°, 1→315°, 2→45°, 3→135°;
        CCW order [0,1,3,4] in vertex-ID space.
        """
        # Apex is the unique degree-4 vertex.
        apex_id = next(
            v for v in range(pyramid_view.num_vertices)
            if pyramid_view.degree(v) == 4
        )
        nbrs = pyramid_view.ordered_neighbors(apex_id)
        # All four base-vertex IDs present, exactly once.
        assert len(nbrs) == 4 and len(set(nbrs)) == 4
        # CCW: each consecutive pair must share a triangular face with apex.
        # Verified by checking the sequence wraps correctly (no repeated IDs).
        # Specific order confirmed from known geometry: [0, 1, 3, 4].
        assert nbrs == [0, 1, 3, 4]

    def test_neighbor_count_equals_degree(self):
        dm = cube()
        view = PlanarEmbeddingView(dm)
        for v in range(view.num_vertices):
            assert len(view.ordered_neighbors(v)) == view.degree(v)

    def test_out_of_range_raises(self):
        view = PlanarEmbeddingView(tetrahedron())
        with pytest.raises(IndexError):
            view.ordered_neighbors(100)


# ---------------------------------------------------------------------------
# PlanarEmbeddingView — degree
# ---------------------------------------------------------------------------


class TestDegree:
    def test_tetrahedron_all_degree_3(self):
        view = PlanarEmbeddingView(tetrahedron())
        for v in range(view.num_vertices):
            assert view.degree(v) == 3

    def test_octahedron_all_degree_4(self):
        view = PlanarEmbeddingView(octahedron())
        for v in range(view.num_vertices):
            assert view.degree(v) == 4

    def test_cube_all_degree_3(self):
        view = PlanarEmbeddingView(cube())
        for v in range(view.num_vertices):
            assert view.degree(v) == 3

    def test_pyramid_apex_degree_4(self, pyramid_view):
        # Apex (label 4, vertex ID 2 in scan order) is the unique degree-4 vertex.
        degrees = [pyramid_view.degree(v) for v in range(pyramid_view.num_vertices)]
        assert degrees.count(4) == 1
        assert max(degrees) == 4

    def test_pyramid_base_degree_3(self, pyramid_view):
        # The four base vertices all have degree 3.
        degrees = [pyramid_view.degree(v) for v in range(pyramid_view.num_vertices)]
        assert degrees.count(3) == 4

    def test_out_of_range_raises(self):
        view = PlanarEmbeddingView(tetrahedron())
        with pytest.raises(IndexError):
            view.degree(-1)


# ---------------------------------------------------------------------------
# PlanarEmbeddingView — face_boundary_vertices
# ---------------------------------------------------------------------------


class TestFaceBoundaryVertices:
    def test_boundary_length_matches_orbit(self):
        dm = cube()
        view = PlanarEmbeddingView(dm)
        face_orbits = dm.face_orbits()
        for fid, orbit in enumerate(face_orbits):
            assert len(view.face_boundary_vertices(fid)) == len(orbit)

    def test_tetrahedron_all_triangles(self):
        view = PlanarEmbeddingView(tetrahedron())
        for fid in range(view.num_faces):
            assert len(view.face_boundary_vertices(fid)) == 3

    def test_boundary_vertices_are_valid_ids(self):
        view = PlanarEmbeddingView(dodecahedron())
        for fid in range(view.num_faces):
            bv = view.face_boundary_vertices(fid)
            assert all(0 <= v < view.num_vertices for v in bv)

    def test_pyramid_outer_face_boundary(self, pyramid_dm, pyramid_view):
        # Outer face has 4 vertices (the base square).
        outer_fid = choose_outer_face(pyramid_dm)
        bv = pyramid_view.face_boundary_vertices(outer_fid)
        assert len(bv) == 4
        # All four boundary vertices are base vertices (degree 3, not the apex).
        for v in bv:
            assert pyramid_view.degree(v) == 3

    def test_out_of_range_raises(self):
        view = PlanarEmbeddingView(tetrahedron())
        with pytest.raises(IndexError):
            view.face_boundary_vertices(100)


# ---------------------------------------------------------------------------
# PlanarEmbeddingView — dart-level access
# ---------------------------------------------------------------------------


class TestDartAccess:
    def test_vertex_of_dart_consistent_with_sigma(self):
        dm = cube()
        view = PlanarEmbeddingView(dm)
        for d in range(dm.num_darts):
            v = view.vertex_of_dart(d)
            # sigma[d] must be at the same vertex
            assert view.vertex_of_dart(dm.sigma[d]) == v

    def test_face_of_dart_consistent_with_phi(self):
        dm = tetrahedron()
        view = PlanarEmbeddingView(dm)
        for d in range(dm.num_darts):
            f = view.face_of_dart(d)
            assert view.face_of_dart(dm.phi(d)) == f

    def test_out_of_range_raises(self):
        view = PlanarEmbeddingView(tetrahedron())
        with pytest.raises(IndexError):
            view.vertex_of_dart(9999)


# ---------------------------------------------------------------------------
# choose_outer_face
# ---------------------------------------------------------------------------


class TestChooseOuterFace:
    @pytest.mark.parametrize(
        "make, expected_boundary_size",
        [
            (tetrahedron, 3),
            (cube, 4),
            (octahedron, 3),
            (dodecahedron, 5),
            (icosahedron, 3),
        ],
    )
    def test_platonic_outer_face_size(self, make, expected_boundary_size):
        dm = make()
        fid = choose_outer_face(dm)
        view = PlanarEmbeddingView(dm)
        assert len(view.face_boundary_vertices(fid)) == expected_boundary_size

    def test_prism4_outer_face_is_quad(self):
        # prism(4) = cube: square caps and square sides; largest face = 4.
        dm = prism(4)
        fid = choose_outer_face(dm)
        view = PlanarEmbeddingView(dm)
        assert len(view.face_boundary_vertices(fid)) == 4

    def test_prism5_outer_face_is_pentagon(self):
        dm = prism(5)
        fid = choose_outer_face(dm)
        view = PlanarEmbeddingView(dm)
        assert len(view.face_boundary_vertices(fid)) == 5

    def test_result_is_valid_face_id(self):
        dm = cube()
        fid = choose_outer_face(dm)
        assert 0 <= fid < len(dm.face_orbits())

    def test_deterministic(self):
        dm = cube()
        assert choose_outer_face(dm) == choose_outer_face(dm)

    def test_pyramid_outer_face_has_4_vertices(self, pyramid_dm):
        fid = choose_outer_face(pyramid_dm)
        view = PlanarEmbeddingView(pyramid_dm)
        assert len(view.face_boundary_vertices(fid)) == 4

    def test_outer_face_is_largest(self):
        # For antiprism(5): triangles (size 3) and pentagons (size 5).
        # Should pick a pentagon.
        dm = antiprism(5)
        fid = choose_outer_face(dm)
        view = PlanarEmbeddingView(dm)
        bv = view.face_boundary_vertices(fid)
        max_size = max(
            len(view.face_boundary_vertices(f))
            for f in range(view.num_faces)
        )
        assert len(bv) == max_size


# ---------------------------------------------------------------------------
# outer_face_anchors
# ---------------------------------------------------------------------------


class TestOuterFaceAnchors:
    def _anchors(self, dm):
        fid = choose_outer_face(dm)
        return outer_face_anchors(dm, fid), fid

    def test_anchors_distinct(self):
        for make in (tetrahedron, cube, octahedron, dodecahedron, icosahedron):
            dm = make()
            (v1, v2, vn), _ = self._anchors(dm)
            assert v1 != v2 and v1 != vn and v2 != vn

    def test_anchors_on_boundary(self):
        dm = cube()
        (v1, v2, vn), fid = self._anchors(dm)
        view = PlanarEmbeddingView(dm)
        bv_set = set(view.face_boundary_vertices(fid))
        assert v1 in bv_set and v2 in bv_set and vn in bv_set

    def test_v1_is_smallest_id(self):
        for make in (tetrahedron, cube, dodecahedron, icosahedron):
            dm = make()
            (v1, v2, vn), fid = self._anchors(dm)
            view = PlanarEmbeddingView(dm)
            bv = view.face_boundary_vertices(fid)
            assert v1 == min(bv)

    def test_v1_v2_adjacent(self):
        dm = cube()
        (v1, v2, vn), _ = self._anchors(dm)
        view = PlanarEmbeddingView(dm)
        assert v2 in view.ordered_neighbors(v1)

    def test_v1_vn_adjacent(self):
        dm = cube()
        (v1, v2, vn), _ = self._anchors(dm)
        view = PlanarEmbeddingView(dm)
        assert vn in view.ordered_neighbors(v1)

    def test_v2_is_predecessor_in_phi_orbit(self):
        """v2 must be the vertex BEFORE v1 in the outer face phi orbit."""
        dm = cube()
        (v1, v2, vn), fid = self._anchors(dm)
        view = PlanarEmbeddingView(dm)
        bv = view.face_boundary_vertices(fid)
        i1 = bv.index(v1)
        assert bv[(i1 - 1) % len(bv)] == v2

    def test_vn_is_successor_in_phi_orbit(self):
        """vn must be the vertex AFTER v1 in the outer face phi orbit."""
        dm = cube()
        (v1, v2, vn), fid = self._anchors(dm)
        view = PlanarEmbeddingView(dm)
        bv = view.face_boundary_vertices(fid)
        i1 = bv.index(v1)
        assert bv[(i1 + 1) % len(bv)] == vn

    def test_triangle_outer_face(self):
        dm = tetrahedron()
        fid = choose_outer_face(dm)
        v1, v2, vn = outer_face_anchors(dm, fid)
        assert v1 != v2 != vn != v1
        view = PlanarEmbeddingView(dm)
        bv_set = set(view.face_boundary_vertices(fid))
        assert {v1, v2, vn} == bv_set

    def test_pyramid_anchors(self, pyramid_dm):
        fid = choose_outer_face(pyramid_dm)
        v1, v2, vn = outer_face_anchors(pyramid_dm, fid)
        view = PlanarEmbeddingView(pyramid_dm)
        bv = view.face_boundary_vertices(fid)
        assert v1 == min(bv)
        assert v2 in view.ordered_neighbors(v1)
        assert vn in view.ordered_neighbors(v1)

    def test_invalid_face_id_raises(self):
        dm = cube()
        with pytest.raises(IndexError):
            outer_face_anchors(dm, 999)

    @pytest.mark.parametrize(
        "make",
        [tetrahedron, cube, octahedron, dodecahedron, icosahedron],
    )
    def test_platonic_anchors_valid(self, make):
        dm = make()
        fid = choose_outer_face(dm)
        v1, v2, vn = outer_face_anchors(dm, fid)
        view = PlanarEmbeddingView(dm)
        bv_set = set(view.face_boundary_vertices(fid))
        assert {v1, v2, vn} <= bv_set
        assert len({v1, v2, vn}) == 3
