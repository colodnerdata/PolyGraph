"""Johnson solid generators.

Parametric generators for pyramid and dipyramid families, plus individual
stubs for Johnson solids that are convex deltahedra but have no simple
parametric form. Stubs raise ``NotImplementedError`` until implemented.

Each function returns ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap

__all__ = [
    "pyramid",
    "dipyramid",
    "bipyramid",
    "snub_disphenoid",
    "triaugmented_triangular_prism",
    "gyroelongated_square_bipyramid",
]


def pyramid(n: int) -> DartMap:
    """Return the n-gonal pyramid.

    Parameters
    ----------
    n : int
        Number of sides in the base polygon.

    Returns
    -------
    DartMap
        Closed combinatorial map with counts
        ``(V, E, F) = (n + 1, 2n, n + 1)``.

    Raises
    ------
    ValueError
        If ``n < 3``.
    """
    if n < 3:
        raise ValueError("pyramid(n) requires n >= 3")

    base = list(range(n))
    apex = n

    faces = [base]
    faces += [[base[(i + 1) % n], base[i], apex] for i in range(n)]
    return DartMap.from_face_lists(
        faces=faces,
        num_vertices=n + 1,
    )


def dipyramid(n: int) -> DartMap:
    """Return the n-gonal dipyramid.

    Parameters
    ----------
    n : int
        Number of vertices in the equatorial cycle.

    Returns
    -------
    DartMap
        Closed combinatorial map with counts ``(V, E, F) = (n + 2, 3n, 2n)``.

    Raises
    ------
    ValueError
        If ``n < 3``.
    """
    if n < 3:
        raise ValueError("dipyramid(n) requires n >= 3")

    equator = list(range(n))
    top = n
    bottom = n + 1

    faces = [[equator[(i + 1) % n], equator[i], top] for i in range(n)]
    faces += [[equator[i], equator[(i + 1) % n], bottom] for i in range(n)]
    return DartMap.from_face_lists(faces=faces, num_vertices=n + 2)


def bipyramid(n: int) -> DartMap:
    """Return the n-gonal bipyramid (alias for :func:`dipyramid`)."""
    return dipyramid(n)


# ---------------------------------------------------------------------------
# Johnson deltahedra without a simple parametric form
# (all faces equilateral triangles; deferred — hardcoded face lists needed)
# ---------------------------------------------------------------------------


def snub_disphenoid() -> DartMap:
    """Return the snub disphenoid (J84).

    Returns
    -------
    DartMap
        Closed combinatorial map with ``(V, E, F) = (8, 18, 12)``.
        All 12 faces are triangles.  Symmetry group D2d.

    Notes
    -----
    4 vertices have degree 5 (0, 2, 4, 6) and 4 have degree 4 (1, 3, 5, 7).
    """
    faces = [
        [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1],
        [1, 5, 6], [1, 6, 2],
        [2, 6, 7], [2, 7, 3],
        [3, 7, 4],
        [4, 7, 6], [4, 6, 5],
    ]
    return DartMap.from_face_lists(faces, num_vertices=8)


def triaugmented_triangular_prism() -> DartMap:
    """Return the triaugmented triangular prism (J51).

    Returns
    -------
    DartMap
        Closed combinatorial map with ``(V, E, F) = (9, 21, 14)``.
        All 14 faces are triangles.  Symmetry group D3h.

    Notes
    -----
    Constructed from a triangular prism (vertices 0–5; top 0,1,2 with
    0↔3, 1↔4, 2↔5) by augmenting each rectangular lateral face with a
    square pyramid: apex 6 over face {0,1,4,3}, apex 7 over {1,2,5,4},
    apex 8 over {2,0,3,5}.
    """
    faces = [
        [0, 1, 2],                          # top cap
        [3, 5, 4],                          # bottom cap
        [0, 6, 1], [1, 6, 4], [4, 6, 3], [3, 6, 0],   # apex 6
        [1, 7, 2], [2, 7, 5], [5, 7, 4], [4, 7, 1],   # apex 7
        [2, 8, 0], [0, 8, 3], [3, 8, 5], [5, 8, 2],   # apex 8
    ]
    return DartMap.from_face_lists(faces, num_vertices=9)


def gyroelongated_square_bipyramid() -> DartMap:
    """Return the gyroelongated square bipyramid (J17).

    Returns
    -------
    DartMap
        Closed combinatorial map with ``(V, E, F) = (10, 24, 16)``.
        All 16 faces are triangles.  Symmetry group D4d.

    Notes
    -----
    Constructed from a square antiprism (vertices 0–7; top ring 0,1,2,3
    CCW from above; bottom ring 4,5,6,7 offset 45°, with bottom[i] between
    top[i] and top[i+1]) with bipyramid apices 8 (top) and 9 (bottom).
    """
    faces = [
        # Square antiprism lateral faces
        [0, 4, 1], [1, 5, 2], [2, 6, 3], [3, 7, 0],   # type-1
        [4, 5, 1], [5, 6, 2], [6, 7, 3], [7, 4, 0],   # type-2
        # Top pyramid (apex 8)
        [0, 1, 8], [1, 2, 8], [2, 3, 8], [3, 0, 8],
        # Bottom pyramid (apex 9)
        [5, 4, 9], [6, 5, 9], [7, 6, 9], [4, 7, 9]
    ]
    return DartMap.from_face_lists(faces, num_vertices=10)
