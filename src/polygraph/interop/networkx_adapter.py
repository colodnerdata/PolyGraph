"""Adapter for converting between DartMap and NetworkX graphs.

This module provides round-trip conversion between PolyGraph's combinatorial
map representation and NetworkX undirected graphs.

The vertex labeling used by :func:`dart_map_to_nx` assigns integer node
indices ``0..V-1`` in the order that vertex orbits are first encountered
during a dart scan from dart 0.  Embedding information (cyclic neighbor
order) is not preserved in the NetworkX graph.

For the reverse direction, :func:`nx_to_dart_map` derives face lists from
the planar half-edge embedding returned by ``networkx.check_planarity``
and then delegates to :meth:`~polygraph.structures.dart_map.DartMap.from_face_lists`.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap

_IMPORT_ERROR_MSG = (
    "networkx is required for graph interoperability. "
    "Install it with: pip install 'polygraph[planar]'"
)


def dart_map_to_nx(dm: DartMap) -> object:
    """Convert a DartMap to an undirected NetworkX graph.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    nx.Graph
        Undirected graph with ``V`` nodes (integers ``0..V-1``) and ``E``
        edges, one per edge orbit of ``dm``.

    Raises
    ------
    ImportError
        If networkx is not installed.
    """
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError(_IMPORT_ERROR_MSG) from exc

    dart_to_vertex: list[int] = [-1] * dm.num_darts
    vertex_id = 0
    for orbit in dm.vertex_orbits():
        for d in orbit:
            dart_to_vertex[d] = vertex_id
        vertex_id += 1

    g: nx.Graph = nx.Graph()
    g.add_nodes_from(range(vertex_id))
    for d in range(dm.num_darts):
        if d < dm.alpha[d]:
            g.add_edge(dart_to_vertex[d], dart_to_vertex[dm.alpha[d]])

    return g


def nx_to_dart_map(g: object) -> DartMap:
    """Convert a connected undirected planar NetworkX graph to a DartMap.

    Uses NetworkX's LR-planarity algorithm to derive a combinatorial planar
    embedding, then extracts all face cycles (including the outer face) via
    half-edge traversal.  Every face of the sphere-like embedding becomes a
    face in the returned DartMap.

    Parameters
    ----------
    g : nx.Graph
        Connected undirected planar graph.  Node labels may be arbitrary
        hashable values; they are mapped to integers ``0..V-1`` internally.

    Returns
    -------
    DartMap
        Combinatorial map whose vertex, edge, and face counts satisfy the
        planar Euler relation ``V - E + F = 2``.

    Raises
    ------
    ImportError
        If networkx is not installed.
    ValueError
        If ``g`` is not planar.
    """
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError(_IMPORT_ERROR_MSG) from exc

    is_planar, embedding = nx.check_planarity(g)
    if not is_planar:
        raise ValueError("Graph is not planar; cannot convert to DartMap.")

    nodes = list(g.nodes())
    node_to_idx: dict = {n: i for i, n in enumerate(nodes)}
    num_vertices = len(nodes)

    # Each directed half-edge belongs to exactly one face in the planar
    # embedding.  Walk each unvisited half-edge to collect its face, marking
    # all half-edges of that face so they are not revisited.
    seen_half_edges: set[tuple] = set()
    faces: list[list[int]] = []
    for u, v in embedding.edges():
        if (u, v) in seen_half_edges:
            continue
        face_nodes = embedding.traverse_face(u, v)
        size = len(face_nodes)
        for i in range(size):
            seen_half_edges.add((face_nodes[i], face_nodes[(i + 1) % size]))
        faces.append([node_to_idx[node] for node in face_nodes])

    return DartMap.from_face_lists(faces, num_vertices)
