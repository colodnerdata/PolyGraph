"""Traversal utilities for combinatorial dart maps.

This module uses a representative-dart convention:
- vertices are represented by a dart in their ``sigma`` orbit,
- faces are represented by a dart in their ``phi`` orbit,
- edges are represented by a dart in their ``alpha`` pair.

For callers that need dense integer identifiers instead of representative
darts, :func:`dart_to_vertex_ids` and :func:`dart_to_face_ids` assign IDs
``0..V-1`` / ``0..F-1`` in orbit scan order.
"""

from __future__ import annotations

from collections.abc import Iterator

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation


def _sigma_permutation(dm: DartMap) -> Permutation:
    """Build the vertex-rotation permutation.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    Permutation
        Permutation for ``sigma``.
    """
    return Permutation.from_sequence(dm.sigma)


def _phi_permutation(dm: DartMap) -> Permutation:
    """Build the face permutation ``phi``.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    Permutation
        Permutation with ``phi[d] = dm.phi(d)``.
    """
    return Permutation.from_sequence([dm.phi(d) for d in range(dm.num_darts)])


def orbit(start: int, perm: Permutation) -> Iterator[int]:
    """Yield one complete orbit under repeated permutation application.

    Parameters
    ----------
    start : int
        Start index in the permutation domain.
    perm : Permutation
        Permutation acting on ``0..n-1``.

    Yields
    ------
    int
        Successive indices in orbit order beginning at ``start`` and ending
        immediately before ``start`` would repeat.

    Raises
    ------
    IndexError
        If ``start`` is outside ``[0, len(perm))``.
    """
    if not 0 <= start < len(perm):
        raise IndexError(
            f"Start index out of range for orbit traversal: {start}."
        )

    i = start
    while True:
        yield i
        i = perm[i]
        if i == start:
            break


def vertex_darts(dm: DartMap, d: int) -> Iterator[int]:
    """Yield darts incident to one vertex in rotational order.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    d : int
        Representative dart of the target vertex orbit.

    Yields
    ------
    int
        Darts in the ``sigma`` orbit of ``d``.

    Raises
    ------
    IndexError
        If ``d`` is outside ``[0, dm.num_darts)``.
    """
    dm.validate_dart(d)
    cur = d
    while True:
        yield cur
        cur = dm.sigma[cur]
        if cur == d:
            break


def face_darts(dm: DartMap, d: int) -> Iterator[int]:
    """Yield darts on one face boundary in face-walk order.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    d : int
        Representative dart of the target face orbit.

    Yields
    ------
    int
        Darts in the ``phi`` orbit of ``d``.

    Raises
    ------
    IndexError
        If ``d`` is outside ``[0, dm.num_darts)``.
    """
    dm.validate_dart(d)
    # Build the face permutation once, then traverse it without repeated
    # DartMap.phi() calls (which perform bounds checks on every step).
    phi_perm = _phi_permutation(dm)
    yield from orbit(d, phi_perm)


def edge_darts(dm: DartMap, d: int) -> Iterator[int]:
    """Yield the two darts of an undirected edge.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    d : int
        Dart index on the edge.

    Yields
    ------
    int
        First ``d``, then ``alpha[d]``.

    Raises
    ------
    IndexError
        If ``d`` is outside ``[0, dm.num_darts)``.
    """
    dm.validate_dart(d)
    yield d
    yield dm.alpha[d]


def all_vertex_orbits(dm: DartMap) -> Iterator[int]:
    """Yield one representative dart for each vertex orbit.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Yields
    ------
    int
        Representative darts in scan order ``0..dm.num_darts-1``.

    Notes
    -----
    Representatives are traversal-derived and not canonical under relabeling.
    """
    sigma = _sigma_permutation(dm)
    seen = [False] * dm.num_darts
    for d in range(dm.num_darts):
        if seen[d]:
            continue
        yield d
        for x in orbit(d, sigma):
            seen[x] = True


def all_face_orbits(dm: DartMap) -> Iterator[int]:
    """Yield one representative dart for each face orbit.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Yields
    ------
    int
        Representative darts in scan order ``0..dm.num_darts-1``.

    Notes
    -----
    Representatives are traversal-derived and not canonical under relabeling.
    """
    phi = _phi_permutation(dm)
    seen = [False] * dm.num_darts
    for d in range(dm.num_darts):
        if seen[d]:
            continue
        yield d
        for x in orbit(d, phi):
            seen[x] = True


def all_edge_orbits(dm: DartMap) -> Iterator[int]:
    """Yield one representative dart for each edge orbit.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Yields
    ------
    int
        Canonical edge representatives satisfying ``d < dm.alpha[d]``.
    """
    for d in range(dm.num_darts):
        if d < dm.alpha[d]:
            yield d


def _scan_order_ids(orbits: list[list[int]], num_darts: int) -> list[int]:
    """Label every dart with the scan-order index of its orbit.

    Parameters
    ----------
    orbits : list[list[int]]
        Disjoint orbits covering ``0..num_darts-1``, in scan order.
    num_darts : int
        Number of darts in the map.

    Returns
    -------
    list[int]
        Array of length ``num_darts`` with ``ids[d] == i`` for every dart
        ``d`` in ``orbits[i]``.
    """
    ids: list[int] = [-1] * num_darts
    for orbit_id, orbit_darts in enumerate(orbits):
        for d in orbit_darts:
            ids[d] = orbit_id
    return ids


def dart_to_vertex_ids(dm: DartMap) -> list[int]:
    """Return a lookup array mapping each dart to its vertex ID.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[int]
        Array of length ``dm.num_darts`` whose entry ``d`` is the vertex ID
        in ``0..V-1`` of the ``sigma`` orbit containing dart ``d``.

    Notes
    -----
    Vertex IDs are assigned in the order that ``sigma`` orbits are first
    encountered when scanning darts from ``0`` upward, i.e. in the order
    returned by ``dm.vertex_orbits()``.  Dart ``0`` therefore always has
    vertex ID ``0``.  This is the convention used by
    :class:`polygraph.algorithms.planar.embedding.PlanarEmbeddingView` and
    by :mod:`polygraph.interop.networkx_adapter`.

    It differs from the *representative-dart* convention used by
    :mod:`polygraph.algorithms.symmetry.orbits`, where a vertex is named by
    the smallest dart in its ``sigma`` orbit rather than by a dense index.

    See Also
    --------
    dart_to_face_ids : Face-ID analogue over ``phi`` orbits.
    all_vertex_orbits : Representative darts of vertex orbits.
    """
    return _scan_order_ids(dm.vertex_orbits(), dm.num_darts)


def dart_to_face_ids(dm: DartMap) -> list[int]:
    """Return a lookup array mapping each dart to its face ID.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[int]
        Array of length ``dm.num_darts`` whose entry ``d`` is the face ID in
        ``0..F-1`` of the ``phi`` orbit containing dart ``d``.

    Notes
    -----
    Face IDs are assigned in the order that ``phi`` orbits are first
    encountered when scanning darts from ``0`` upward, i.e. in the order
    returned by ``dm.face_orbits()``.  Dart ``0`` therefore always has face
    ID ``0``.  This is the convention used by
    :class:`polygraph.algorithms.planar.embedding.PlanarEmbeddingView` and
    by :mod:`polygraph.interop.networkx_adapter`.

    It differs from the *representative-dart* convention used by
    :mod:`polygraph.algorithms.symmetry.orbits`, where a face is named by
    the smallest dart in its ``phi`` orbit rather than by a dense index.

    See Also
    --------
    dart_to_vertex_ids : Vertex-ID analogue over ``sigma`` orbits.
    all_face_orbits : Representative darts of face orbits.
    """
    return _scan_order_ids(dm.face_orbits(), dm.num_darts)


def neighbors(dm: DartMap, vertex: int) -> Iterator[int]:
    """Yield adjacent vertex representatives around a vertex.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    vertex : int
        Representative dart of the source vertex.

    Yields
    ------
    int
        One representative dart in each adjacent vertex, in the cyclic order
        induced by ``vertex_darts(dm, vertex)``.

    Notes
    -----
    Multiplicity is preserved: parallel edges emit repeated neighbors.
    """
    for d in vertex_darts(dm, vertex):
        yield dm.alpha[d]


def vertices_of_face(dm: DartMap, face: int) -> Iterator[int]:
    """Yield boundary vertex representatives of a face.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    face : int
        Representative dart of the face.

    Yields
    ------
    int
        Boundary vertices represented by darts, in face boundary order.

    Notes
    -----
    Under the representative-dart model, this is the same sequence as
    ``face_darts(dm, face)``.
    """
    yield from face_darts(dm, face)


def faces_incident_to_vertex(dm: DartMap, vertex: int) -> Iterator[int]:
    """Yield incident face representatives around a vertex.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    vertex : int
        Representative dart of the vertex.

    Yields
    ------
    int
        Incident faces represented by darts, in rotational order.

    Notes
    -----
    Under the representative-dart model, this is the same sequence as
    ``vertex_darts(dm, vertex)``.
    """
    yield from vertex_darts(dm, vertex)


def adjacent_face_pairs(dm: DartMap) -> Iterator[tuple[int, int]]:
    """Yield adjacent face representative pairs across edges.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Yields
    ------
    tuple[int, int]
        Pairs ``(face_a, face_b)`` as ``(e, dm.alpha[e])`` for each
        representative edge dart ``e``.

    Notes
    -----
    Exactly one pair is yielded per undirected edge.
    """
    for e in all_edge_orbits(dm):
        yield (e, dm.alpha[e])
