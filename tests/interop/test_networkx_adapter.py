"""Tests for the NetworkX interop adapter."""

from __future__ import annotations

import pytest

from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.interop.networkx_adapter import dart_map_to_nx, nx_to_dart_map

nx = pytest.importorskip("networkx")

# (generator, expected_V, expected_E, expected_F)
PLATONIC = [
    (tetrahedron, 4, 6, 4),
    (cube, 8, 12, 6),
    (octahedron, 6, 12, 8),
    (dodecahedron, 20, 30, 12),
    (icosahedron, 12, 30, 20),
]


class TestDartMapToNx:
    def test_cube_node_and_edge_counts(self):
        g = dart_map_to_nx(cube())
        assert g.number_of_nodes() == 8
        assert g.number_of_edges() == 12

    @pytest.mark.parametrize("make, v, e, f", PLATONIC)
    def test_platonic_node_and_edge_counts(self, make, v, e, f):
        g = dart_map_to_nx(make())
        assert g.number_of_nodes() == v
        assert g.number_of_edges() == e

    def test_nodes_are_integers(self):
        g = dart_map_to_nx(cube())
        assert all(isinstance(n, int) for n in g.nodes())

    def test_nodes_are_zero_indexed(self):
        g = dart_map_to_nx(cube())
        assert set(g.nodes()) == set(range(8))

    def test_result_is_networkx_graph(self):
        g = dart_map_to_nx(tetrahedron())
        assert isinstance(g, nx.Graph)

    def test_result_is_planar(self):
        for make, *_ in PLATONIC:
            g = dart_map_to_nx(make())
            assert nx.check_planarity(g)[0]


class TestNxToDartMap:
    @pytest.mark.parametrize("make, v, e, f", PLATONIC)
    def test_round_trip_vef_counts(self, make, v, e, f):
        dm = make()
        g = dart_map_to_nx(dm)
        dm2 = nx_to_dart_map(g)

        assert len(dm2.vertex_orbits()) == v
        assert dm2.num_edges == e
        assert len(dm2.face_orbits()) == f

    @pytest.mark.parametrize("make, v, e, f", PLATONIC)
    def test_euler_characteristic(self, make, v, e, f):
        dm = nx_to_dart_map(dart_map_to_nx(make()))
        assert dm.euler_characteristic() == 2

    def test_non_planar_raises(self):
        k5 = nx.complete_graph(5)
        with pytest.raises(ValueError, match="not planar"):
            nx_to_dart_map(k5)

    def test_arbitrary_node_labels(self):
        g = nx.Graph()
        g.add_edges_from([("a", "b"), ("b", "c"), ("c", "a")])
        dm = nx_to_dart_map(g)
        assert len(dm.vertex_orbits()) == 3
        assert dm.num_edges == 3
        assert len(dm.face_orbits()) == 2
