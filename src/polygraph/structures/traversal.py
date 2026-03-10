"""Traversal utilities for combinatorial dart maps.

This module uses a representative-dart convention:
- vertices are represented by a dart in their ``sigma`` orbit,
- faces are represented by a dart in their ``phi`` orbit,
- edges are represented by a dart in their ``alpha`` pair.
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
        raise IndexError(f"Start index out of range for orbit traversal: {start}.")

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
    if not 0 <= d < dm.num_darts:
        raise IndexError(
            f"Dart index {d} out of range [0, {dm.num_darts})."
        )
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
    if not 0 <= d < dm.num_darts:
        raise IndexError(
            f"Dart index {d} out of range [0, {dm.num_darts})."
        )
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
    if not 0 <= d < dm.num_darts:
        raise IndexError(
            f"Dart index {d} out of range [0, {dm.num_darts})."
        )
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
