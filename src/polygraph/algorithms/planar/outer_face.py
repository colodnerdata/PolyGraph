"""Outer face selection and anchor vertex computation for planar drawing.

These functions implement the two setup steps required before any canonical
ordering algorithm can run:

1. :func:`choose_outer_face` — select the face to treat as the unbounded
   (outer) face.  Heuristic: largest boundary; deterministic tie-break.

2. :func:`outer_face_anchors` — identify the three anchor vertices
   ``(v1, v2, vn)`` on the chosen outer face, following the conventions of
   de Fraysseix, Pach & Pollack (1990) and Kant (1996):

   - ``v1`` is the outer-face vertex with the smallest integer ID.
   - ``v2`` is the **CCW** neighbour of ``v1`` on the outer boundary, i.e.,
     the vertex immediately **before** ``v1`` in the phi orbit (which is CW
     in the plane for the outer face).
   - ``vn`` is the **CW** neighbour of ``v1`` on the outer boundary, i.e.,
     the vertex immediately **after** ``v1`` in the phi orbit.
   - Both ``(v1, v2)`` and ``(v1, vn)`` are edges of the outer face.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.traversal import face_darts


def _build_dart_to_vertex(dm: DartMap) -> list[int]:
    """Return a dart → vertex-ID lookup array.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[int]
        Array of length ``dm.num_darts``.
    """
    dart_to_vertex: list[int] = [-1] * dm.num_darts
    for vid, orbit in enumerate(dm.vertex_orbits()):
        for d in orbit:
            dart_to_vertex[d] = vid
    return dart_to_vertex


def choose_outer_face(dm: DartMap) -> int:
    """Select the outer face for planar drawing.

    Chooses the face with the most boundary vertices.  Among ties, the face
    with the lowest face ID (equivalently, the smallest representative dart)
    is returned.  This is deterministic for all inputs including Platonic
    solids where every face has the same size.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    int
        Face ID in ``[0, F)`` where ``F`` is the number of faces.
    """
    face_orbits = dm.face_orbits()
    best_fid = 0
    best_size = len(face_orbits[0])
    for fid in range(1, len(face_orbits)):
        size = len(face_orbits[fid])
        if size > best_size:
            best_size = size
            best_fid = fid
    return best_fid


def outer_face_anchors(
    dm: DartMap, outer_face: int
) -> tuple[int, int, int]:
    """Return the three anchor vertices for canonical-ordering algorithms.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    outer_face : int
        Face ID of the chosen outer face, as returned by
        :func:`choose_outer_face`.

    Returns
    -------
    v1 : int
        Outer-face vertex with the smallest vertex ID.
    v2 : int
        CCW neighbour of ``v1`` on the outer boundary (predecessor in phi
        orbit = CCW in the plane drawing).
    vn : int
        CW neighbour of ``v1`` on the outer boundary (successor in phi
        orbit = CW in the plane drawing).

    Raises
    ------
    IndexError
        If ``outer_face`` is out of range.
    ValueError
        If the outer face has fewer than 3 vertices.

    Notes
    -----
    Both ``(v1, v2)`` and ``(v1, vn)`` are edges of the outer face.  All
    three vertices are distinct (guaranteed for any face of size ≥ 3).
    """
    face_orbits = dm.face_orbits()
    if not 0 <= outer_face < len(face_orbits):
        raise IndexError(
            f"Face ID {outer_face} out of range [0, {len(face_orbits)})."
        )

    orbit = face_orbits[outer_face]
    k = len(orbit)
    if k < 3:
        raise ValueError(
            f"Outer face has {k} vertices; at least 3 are required."
        )

    dart_to_vertex = _build_dart_to_vertex(dm)

    # Boundary vertices in phi orbit order (CW in the plane for outer face).
    bv = [dart_to_vertex[d] for d in face_darts(dm, orbit[0])]

    # v1 = vertex with smallest ID on the boundary.
    min_vid = min(bv)
    i1 = bv.index(min_vid)

    v1 = bv[i1]
    v2 = bv[(i1 - 1) % k]   # predecessor in phi = CCW neighbour in plane
    vn = bv[(i1 + 1) % k]   # successor in phi   = CW  neighbour in plane

    return v1, v2, vn
