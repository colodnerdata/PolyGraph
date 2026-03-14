"""Generators for selected Platonic solids via explicit face cycles.

The generators in this module construct combinatorial maps directly from
oriented face-vertex cycles using :meth:`polygraph.structures.dart_map.DartMap.
from_face_lists`.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap

__all__ = ["tetrahedron", "dodecahedron", "icosahedron"]


def tetrahedron() -> DartMap:
    """Build the tetrahedron map.

    Returns
    -------
    DartMap
        Tetrahedron map with counts ``(V, E, F) = (4, 6, 4)``.

    Notes
    -----
    The 4 triangular faces are oriented consistently so each undirected edge
    appears exactly twice with opposite directions.
    """
    faces = [
        [0, 1, 2],
        [0, 3, 1],
        [1, 3, 2],
        [0, 2, 3],
    ]
    return DartMap.from_face_lists(faces=faces, num_vertices=4)


def dodecahedron() -> DartMap:
    """Build the dodecahedron map.

    Returns
    -------
    DartMap
        Dodecahedron map with counts ``(V, E, F) = (20, 30, 12)``.

    Notes
    -----
    The 12 pentagonal faces are oriented consistently so each undirected edge
    appears exactly twice with opposite directions.
    """
    faces = [
        [0, 4, 3, 2, 1],
        [1, 2, 9, 19, 5],
        [7, 16, 11, 12, 17],
        [10, 14, 13, 12, 11],
        [6, 15, 10, 11, 16],
        [0, 1, 5, 15, 6],
        [8, 17, 12, 13, 18],
        [2, 3, 8, 18, 9],
        [9, 18, 13, 14, 19],
        [5, 19, 14, 10, 15],
        [3, 4, 7, 17, 8],
        [0, 6, 16, 7, 4],
    ]
    return DartMap.from_face_lists(faces=faces, num_vertices=20)


def icosahedron() -> DartMap:
    """Build the icosahedron map.

    Returns
    -------
    DartMap
        Icosahedron map with counts ``(V, E, F) = (12, 30, 20)``.

    Notes
    -----
    The 20 triangular faces are oriented consistently so each undirected edge
    appears exactly twice with opposite directions.
    """
    faces = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ]
    return DartMap.from_face_lists(faces=faces, num_vertices=12)
