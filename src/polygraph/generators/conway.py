"""Conway operators on combinatorial maps, with provenance.

Every operator takes a :class:`~polygraph.structures.dart_map.DartMap` and
returns a :class:`ConwayResult`: the new map together with *provenance*, which
records the old object (vertex, edge, face or dart) each new vertex, dart and
face was built from, and which raw vertex keys were identified into each new
vertex.

Concepts
--------
A ``DartMap`` is a pair of permutations ``(sigma, alpha)`` on darts
(oriented half-edges).  ``sigma(d)`` is the next dart counterclockwise around
the tail of ``d``; ``alpha(d)`` is the opposite dart of the same edge; and
``phi = sigma^{-1} o alpha`` is the next dart along the face to the *left* of
``d``.  Vertices are the cycles of ``sigma``, faces the cycles of ``phi``, and
edges the pairs ``{d, alpha(d)}``.

A Conway operator is described geometrically ("put a vertex at every edge
midpoint").  Translated to darts, such a description produces *raw* vertices,
and the same geometric vertex may be reached from several darts, edges or
faces; those raw vertices must be *identified* (merged).  The ``merged_from``
field of :class:`ConwayResult` records exactly which raw keys were merged.
Canonical ids of old objects are the smallest dart of their orbit, so they do
not depend on how orbits happen to be listed.

Orientation convention
----------------------
Faces are assumed counterclockwise when viewed from outside, so ``sigma``
turns counterclockwise about a vertex.  Chirality only matters for operators
that are not mirror-symmetric (``snub``, ``whirl``); :func:`dual` is
achiral.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from polygraph.structures.dart_map import DartMap
from polygraph.structures.dual import dual_of

__all__ = ["ConwayResult", "Origin", "dual"]

OriginKind = Literal["vertex", "edge", "face", "dart"]
"""The kind of old object a new object was built from."""

RawKey = tuple[str, int]
"""A raw vertex key such as ``("face", 4)``: (kind, canonical old id)."""


@dataclass(frozen=True)
class Origin:
    """The old object a new vertex, dart or face was built from.

    Attributes
    ----------
    kind : {"vertex", "edge", "face", "dart"}
        Kind of the old object.
    old_id : int
        Canonical id of the old object: the dart itself for ``"dart"``, and
        the smallest dart of its orbit for ``"vertex"`` and ``"face"``
        (``min(d, alpha(d))`` for ``"edge"``).
    role : str
        What the new object is relative to the old one (for example
        ``"dual"``), so an operator that builds several new objects from one
        old object can tell them apart.
    """

    kind: OriginKind
    old_id: int
    role: str = ""


@dataclass(frozen=True)
class ConwayResult:
    """The output of a Conway operator.

    All sequences are indexed like the corresponding orbit lists of ``dm``:
    ``vertex_origin[i]`` describes ``dm.vertex_orbits()[i]``,
    ``face_origin[j]`` describes ``dm.face_orbits()[j]``, and
    ``dart_origin[d]`` describes dart ``d``.

    Attributes
    ----------
    dm : DartMap
        The new map.
    vertex_origin : tuple[Origin, ...]
        Old object each new vertex was built from.
    dart_origin : tuple[Origin | None, ...]
        Old object each new dart was built from; ``None`` for a dart with no
        old counterpart.
    face_origin : tuple[Origin, ...]
        Old object each new face was built from.
    merged_from : tuple[tuple[RawKey, ...], ...]
        For each new vertex, the raw keys identified into it.  A length-1
        tuple means nothing was merged.
    """

    dm: DartMap
    vertex_origin: tuple[Origin, ...]
    dart_origin: tuple[Origin | None, ...]
    face_origin: tuple[Origin, ...]
    merged_from: tuple[tuple[RawKey, ...], ...]


def _canonical_ids(orbits: list[list[int]], num_darts: int) -> list[int]:
    """Map every dart to the smallest dart of its orbit.

    Parameters
    ----------
    orbits : list[list[int]]
        Disjoint orbits covering all darts ``0..num_darts-1``.
    num_darts : int
        Total number of darts.

    Returns
    -------
    list[int]
        ``ids[d] == min(orbit containing d)``.
    """
    ids = [0] * num_darts
    for orbit in orbits:
        smallest = min(orbit)
        for d in orbit:
            ids[d] = smallest
    return ids


def dual(dm: DartMap) -> ConwayResult:
    """Return the dual of a map: faces become vertices and vice versa.

    Parameters
    ----------
    dm : DartMap
        Any closed oriented map (in particular the map of a convex genus-0
        polyhedron).

    Returns
    -------
    ConwayResult
        ``result.dm`` is the dual map ``(phi, alpha)`` on the *same* darts.
        Provenance:

        - ``vertex_origin``: one ``Origin("face", f, "dual")`` per new
          vertex, where ``f`` is the smallest dart of the old face.
        - ``face_origin``: one ``Origin("vertex", v, "dual")`` per new face,
          where ``v`` is the smallest dart of the old vertex.
        - ``dart_origin[d] == Origin("dart", d, "dual")``.
        - ``merged_from``: singletons ``(("face", f),)``; the dual merges
          nothing.

    Notes
    -----
    Let ``phi = sigma^{-1} o alpha`` be the face permutation.  The dual keeps
    the darts and ``alpha`` and replaces ``sigma`` by ``phi``:

    ``(sigma*, alpha*) = (phi, alpha)``.

    **New vertices.**  The vertex orbits of the dual are the cycles of
    ``sigma* = phi``, which are exactly the old faces.  Dart ``d`` therefore
    *starts at* the old face ``face(d)`` to its left, and the cyclic order of
    darts leaving a new vertex is the counterclockwise boundary order of the
    old face.

    **New faces.**  The face permutation of the dual is

    ``phi* = (sigma*)^{-1} o alpha = phi^{-1} o alpha``
    ``     = alpha o sigma o alpha``,

    using ``phi^{-1} = alpha^{-1} o sigma = alpha o sigma``.  So ``phi*`` is
    ``sigma`` conjugated by ``alpha``: its cycles are the sets
    ``alpha(orbit)`` of old vertex orbits.  The new face for an old vertex
    ``v`` is the set of darts *pointing into* ``v``, equivalently the dart
    ``d`` has old vertex ``head(d)`` on its left.

    **Edges.**  An old edge ``{d, alpha(d)}`` becomes the new edge with the
    same two darts; it crosses the old edge, joining ``face(d)`` to
    ``face(alpha(d))``.

    **Orientation.**  Around the old face ``f`` the darts leaving the new
    vertex ``f`` appear in counterclockwise order, so ``sigma*`` still turns
    counterclockwise and the dual keeps the outward-counterclockwise
    convention.  Hence the dual is orientation-preserving and achiral:
    ``dual(mirror(M))`` is ``mirror(dual(M))`` up to relabelling by
    ``alpha``.

    **Counts.**  ``(V', E', F') = (F, E, V)``; vertex degrees and face sizes
    swap.  Convexity is preserved: the dual of a 3-connected planar graph is
    3-connected planar (polarity).

    **Involution up to relabelling.**  ``dual(dual(M))`` equals ``M``
    conjugated by ``alpha``, *not* ``M`` itself: its ``sigma`` is
    ``alpha o sigma o alpha``.  The correct check is
    ``is_isomorphism(M, dual(dual(M)).dm, M.alpha)``.

    Complexity is ``O(n)`` in the number of darts.

    Examples
    --------
    >>> from polygraph.generators.platonic import cube
    >>> result = dual(cube())
    >>> len(result.dm.vertex_orbits()), len(result.dm.face_orbits())
    (6, 8)
    """
    result = dual_of(dm)
    old_vertex_id = _canonical_ids(dm.vertex_orbits(), dm.num_darts)

    # Dual vertex orbits are old face orbits (sigma* = phi), so min(orbit) is
    # the old face's canonical id whatever order the orbits are listed in.
    vertex_ids = [min(orbit) for orbit in result.vertex_orbits()]
    vertex_origin = tuple(Origin("face", f, "dual") for f in vertex_ids)
    merged_from = tuple((("face", f),) for f in vertex_ids)

    # Every dart of a dual face points into the same old vertex, so any one
    # dart identifies it: orbit[0] points into the old vertex containing
    # alpha(orbit[0]).
    face_origin = tuple(
        Origin("vertex", old_vertex_id[dm.alpha[orbit[0]]], "dual")
        for orbit in result.face_orbits()
    )

    dart_origin = tuple(Origin("dart", d, "dual") for d in range(dm.num_darts))
    return ConwayResult(
        dm=result,
        vertex_origin=vertex_origin,
        dart_origin=dart_origin,
        face_origin=face_origin,
        merged_from=merged_from,
    )
