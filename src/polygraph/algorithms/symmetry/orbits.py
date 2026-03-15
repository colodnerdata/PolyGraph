"""Orbit decomposition under dart-map automorphism groups.

Given a generating set for ``Aut(dm)``, the functions in this module
partition darts, vertices, edges, and faces into equivalence classes under
the group action.

All functions use **union-find with path compression** and require only a
single pass over the generators — the union-find closure is equivalent to
the full group orbit partition.
"""

from __future__ import annotations

from collections import defaultdict

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation

# ---------------------------------------------------------------------------
# Internal union-find helpers
# ---------------------------------------------------------------------------


def _make_uf(n: int) -> list[int]:
    """Return a fresh union-find parent array of size *n*."""
    return list(range(n))


def _find(parent: list[int], x: int) -> int:
    """Return the representative of *x* with path-halving compression."""
    while parent[x] != x:
        parent[x] = parent[parent[x]]  # path halving
        x = parent[x]
    return x


def _union(parent: list[int], a: int, b: int) -> None:
    """Merge the sets containing *a* and *b* in place."""
    ra, rb = _find(parent, a), _find(parent, b)
    if ra != rb:
        parent[ra] = rb


def _collect_orbits(parent: list[int]) -> list[list[int]]:
    """Convert a union-find parent array into a list of orbit lists."""
    groups: dict[int, list[int]] = defaultdict(list)
    for i in range(len(parent)):
        groups[_find(parent, i)].append(i)
    return list(groups.values())


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_orbits(
    generators: list[Permutation], n: int
) -> list[list[int]]:
    """Partition ``0..n-1`` into orbits under the generated group.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for the symmetry group.
    n : int
        Domain size.

    Returns
    -------
    list[list[int]]
        Disjoint lists whose union is ``range(n)``.
    """
    parent = _make_uf(n)
    for g in generators:
        for i in range(n):
            _union(parent, i, g[i])
    return _collect_orbits(parent)


def dart_orbits(
    generators: list[Permutation], dm: DartMap
) -> list[list[int]]:
    """Return orbits of dart indices.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for ``Aut(dm)``.
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[list[int]]
        Each inner list is one dart orbit.
    """
    return compute_orbits(generators, dm.num_darts)


def vertex_orbits(
    generators: list[Permutation], dm: DartMap
) -> list[list[int]]:
    """Return orbits of vertices under the automorphism group.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for ``Aut(dm)``.
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[list[int]]
        Each inner list contains the representative darts of one vertex orbit.

    Notes
    -----
    Vertices are represented by the smallest dart in their ``sigma`` orbit
    (the representative-dart convention from
    :mod:`polygraph.structures.traversal`).
    """
    # Map each dart to its vertex representative (min dart in sigma orbit).
    dart_to_vertex: list[int] = [0] * dm.num_darts
    for cycle in dm.vertex_orbits():
        rep = min(cycle)
        for d in cycle:
            dart_to_vertex[d] = rep

    # Collect distinct vertex representatives.
    vertex_reps = sorted(set(dart_to_vertex))
    rep_to_idx = {r: i for i, r in enumerate(vertex_reps)}
    num_vertices = len(vertex_reps)

    parent = _make_uf(num_vertices)
    for g in generators:
        for d in range(dm.num_darts):
            v_src = rep_to_idx[dart_to_vertex[d]]
            v_dst = rep_to_idx[dart_to_vertex[g[d]]]
            _union(parent, v_src, v_dst)

    raw = _collect_orbits(parent)
    # Convert index orbits back to representative darts.
    return [[vertex_reps[i] for i in orbit] for orbit in raw]


def edge_orbits(
    generators: list[Permutation], dm: DartMap
) -> list[list[int]]:
    """Return orbits of undirected edges under the automorphism group.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for ``Aut(dm)``.
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[list[int]]
        Each inner list contains the canonical darts of one edge orbit.
        The canonical dart of an edge is ``min(d, alpha[d])``.
    """
    alpha = dm.alpha
    # Canonical dart for each dart's edge.
    dart_to_edge: list[int] = [min(d, alpha[d]) for d in range(dm.num_darts)]

    edge_reps = sorted(d for d in range(dm.num_darts) if d < alpha[d])
    rep_to_idx = {r: i for i, r in enumerate(edge_reps)}
    num_edges = len(edge_reps)

    parent = _make_uf(num_edges)
    for g in generators:
        for d in range(dm.num_darts):
            if d < alpha[d]:  # process each edge once
                e_src = rep_to_idx[dart_to_edge[d]]
                # Image of edge (d, alpha[d]) under g.
                e_dst = rep_to_idx[min(g[d], g[alpha[d]])]
                _union(parent, e_src, e_dst)

    raw = _collect_orbits(parent)
    return [[edge_reps[i] for i in orbit] for orbit in raw]


def face_orbits(
    generators: list[Permutation], dm: DartMap
) -> list[list[int]]:
    """Return orbits of faces under the automorphism group.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for ``Aut(dm)``.
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[list[int]]
        Each inner list contains the representative darts of one face orbit.

    Notes
    -----
    Faces are represented by the smallest dart in their ``phi`` orbit.
    """
    # Map each dart to its face representative (smallest dart in phi orbit).
    dart_to_face: list[int] = [0] * dm.num_darts
    for cycle in dm.face_orbits():
        rep = min(cycle)
        for d in cycle:
            dart_to_face[d] = rep

    face_reps = sorted(set(dart_to_face))
    rep_to_idx = {r: i for i, r in enumerate(face_reps)}
    num_faces = len(face_reps)

    # For orientation-reversing generators, the image of a face orbit may be
    # traversed in reverse (phi⁻¹ instead of phi).  Either way, all darts of
    # the image face share the same face representative, so we just look up
    # dart_to_face[g[d]] for any dart d in the face.
    parent = _make_uf(num_faces)
    for g in generators:
        for d in range(dm.num_darts):
            f_src = rep_to_idx[dart_to_face[d]]
            f_dst = rep_to_idx[dart_to_face[g[d]]]
            _union(parent, f_src, f_dst)

    raw = _collect_orbits(parent)
    return [[face_reps[i] for i in orbit] for orbit in raw]
