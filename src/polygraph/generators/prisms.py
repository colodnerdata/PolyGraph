"""Prism and antiprism generators.

Parametric generators for the regular prism and antiprism families. Vertex
labels alternate between the top and bottom caps: even indices belong to the
top cycle and odd indices belong to the bottom cycle.

Each function returns ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap


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
    """
    if n < 3:
        raise ValueError("prism(n) requires n >= 3")

    top = list(range(0, 2 * n, 2))
    bottom = list(range(1, 2 * n, 2))
    faces = [top, bottom[::-1]]
    # The top cap contributes each top edge as top[i] -> top[i + 1], while the
    # side quad uses the same undirected edge in reverse as
    # top[i + 1] -> top[i]. The bottom cap is reversed for the same reason, so
    # every cap edge appears exactly twice with opposite directions.
    faces += [
        [top[i], bottom[i], bottom[(i + 1) % n], top[(i + 1) % n]]
        for i in range(n)
    ]
    # Each vertical edge appears in two adjacent side quads:
    # top[i] -> bottom[i] in quad i and bottom[i] -> top[i] in quad i - 1.
    # Each bottom edge appears as bottom[i] -> bottom[i + 1] here and
    # bottom[i + 1] -> bottom[i] in the reversed bottom cap.
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
    """
    if n < 3:
        raise ValueError("antiprism(n) requires n >= 3")

    top = list(range(0, 2 * n, 2))
    bottom = list(range(1, 2 * n, 2))
    faces = [top, bottom[::-1]]
    # The first triangle family uses each top edge as
    # top[i + 1] -> top[i], opposite the top-cap direction
    # top[i] -> top[i + 1]. Its slanted edges are then paired by the adjacent
    # triangles in the second family.
    faces += [
        [top[i], bottom[i], top[(i + 1) % n]]
        for i in range(n)
    ]
    # The second triangle family completes the band so that
    # top[i + 1] -> bottom[i] reverses bottom[i] -> top[i + 1] from the first
    # family, and bottom[i] -> bottom[i + 1] reverses the direction used in the
    # reversed bottom cap. This makes every undirected edge appear exactly
    # twice, once per adjacent face, with opposite winding.
    faces += [
        [top[(i + 1) % n], bottom[i], bottom[(i + 1) % n]]
        for i in range(n)
    ]
    return DartMap.from_face_lists(faces, num_vertices=2 * n)
