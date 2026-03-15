"""Prism and antiprism generators.

Parametric generators for the regular prism and antiprism families. Vertex
labels alternate between the top and bottom caps: even indices belong to the
top cycle and odd indices belong to the bottom cycle.

Each function returns ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap

__all__ = ["prism", "antiprism"]


def prism(n: int) -> DartMap:
    """Return the n-gonal prism.

    Parameters
    ----------
    n : int
        Number of sides in each congruent polygonal cap.

    Returns
    -------
    DartMap
        Closed combinatorial map of the prism with ``2 * n`` vertices.

    Raises
    ------
    ValueError
        If ``n < 3``.

    Notes
    -----
    The top cap traverses each top edge as ``top[i] -> top[i + 1]``, and the
    corresponding side quad traverses the same undirected edge in reverse as
    ``top[i + 1] -> top[i]``. The bottom cap is reversed for the same reason.

    Each vertical edge appears in two adjacent side quads, once as
    ``top[i] -> bottom[i]`` and once as ``bottom[i] -> top[i]``. Each bottom
    edge appears once in a side quad and once in the reversed bottom cap with
    opposite directions, so every undirected edge appears exactly twice with
    opposite winding across adjacent faces.
    """
    if n < 3:
        raise ValueError("prism(n) requires n >= 3")

    top = list(range(0, 2 * n, 2))
    bottom = list(range(1, 2 * n, 2))
    faces = [top, bottom[::-1]]
    faces += [
        [top[i], bottom[i], bottom[(i + 1) % n], top[(i + 1) % n]]
        for i in range(n)
    ]
    return DartMap.from_face_lists(faces, num_vertices=2 * n)


def antiprism(n: int) -> DartMap:
    """Return the n-gonal antiprism.

    Parameters
    ----------
    n : int
        Number of sides in each congruent polygonal cap.

    Returns
    -------
    DartMap
        Closed combinatorial map of the antiprism with ``2 * n`` vertices.

    Raises
    ------
    ValueError
        If ``n < 3``.

    Notes
    -----
    The first triangle family traverses each top edge as
    ``top[i + 1] -> top[i]``, opposite the top-cap direction
    ``top[i] -> top[i + 1]``. Its slanted edges are paired by adjacent
    triangles in the second family.

    The second triangle family completes the band so that
    ``top[i + 1] -> bottom[i]`` reverses ``bottom[i] -> top[i + 1]`` from the
    first family, while ``bottom[i] -> bottom[i + 1]`` reverses the direction
    used in the reversed bottom cap. Every undirected edge therefore appears
    exactly twice, once per adjacent face, with opposite winding.
    """
    if n < 3:
        raise ValueError("antiprism(n) requires n >= 3")

    top = list(range(0, 2 * n, 2))
    bottom = list(range(1, 2 * n, 2))
    faces = [top, bottom[::-1]]
    faces += [
        [top[i], bottom[i], top[(i + 1) % n]]
        for i in range(n)
    ]
    faces += [
        [top[(i + 1) % n], bottom[i], bottom[(i + 1) % n]]
        for i in range(n)
    ]
    return DartMap.from_face_lists(faces, num_vertices=2 * n)
