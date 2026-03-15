"""Symmetry analysis for combinatorial dart maps.

This package computes the automorphism group of a
:class:`~polygraph.structures.dart_map.DartMap` and partitions its
structural elements (darts, vertices, edges, faces) into symmetry orbits.

Automorphism computation delegates to pynauty (which uses nauty internally).
Install the optional dependency with ``pip install 'polygraph[symmetry]'``.
"""

from polygraph.algorithms.symmetry.automorphisms import (
    automorphism_group_order,
    compute_automorphism_generators,
    is_orientation_preserving,
)
from polygraph.algorithms.symmetry.orbits import (
    compute_orbits,
    dart_orbits,
    edge_orbits,
    face_orbits,
    vertex_orbits,
)

__all__ = [
    "automorphism_group_order",
    "compute_automorphism_generators",
    "is_orientation_preserving",
    "compute_orbits",
    "dart_orbits",
    "vertex_orbits",
    "edge_orbits",
    "face_orbits",
]
