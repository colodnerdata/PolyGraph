"""Prism and antiprism generators.

Parametric generators for the regular prism and antiprism families. Vertex
labels alternate between the top and bottom caps: even indices belong to the
top cycle and odd indices belong to the bottom cycle.

Each function returns ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap


def prism(n: int) -> DartMap:
    """Return the n-gonal prism."""
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
    """Return the n-gonal antiprism."""
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
