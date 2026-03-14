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
        Closed combinatorial map with counts ``(V, E, F) = (n + 1, 2n, n + 1)``.

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
    return DartMap.from_face_lists(faces=faces, num_vertices=n + 1)


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


def snub_disphenoid():
    """Return the snub disphenoid (J84).  V=8, E=18, F=12."""
    raise NotImplementedError


def triaugmented_triangular_prism():
    """Return the triaugmented triangular prism (J51).  V=9, E=21, F=14."""
    raise NotImplementedError


def gyroelongated_square_bipyramid():
    """Return the gyroelongated square bipyramid (J17).  V=10, E=24, F=16."""
    raise NotImplementedError
