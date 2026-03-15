"""Symmetry analysis for combinatorial dart maps.

This package computes the automorphism group of a
:class:`~polygraph.structures.dart_map.DartMap`, partitions its structural
elements (darts, vertices, edges, faces) into symmetry orbits, and classifies
the automorphism group as a named point group.

Automorphism computation delegates to pynauty (which uses nauty internally).
Install the optional dependency with ``pip install 'polygraph[symmetry]'``.
"""
from __future__ import annotations

from polygraph.algorithms.symmetry.automorphisms import (
    automorphism_group_order,
    compute_automorphism_generators,
    is_orientation_preserving,
)
from polygraph.algorithms.symmetry.classify import (
    UnknownSymmetryError,
    classify_symmetry,
)
from polygraph.algorithms.symmetry.orbits import (
    compute_orbits,
    dart_orbits,
    edge_orbits,
    face_orbits,
    vertex_orbits,
)
from polygraph.algorithms.symmetry.point_groups import (
    PointGroup,
    cyclic,
    dihedral,
    icosahedral,
    octahedral,
    tetrahedral,
)

__all__ = [
    # automorphisms
    "automorphism_group_order",
    "compute_automorphism_generators",
    "is_orientation_preserving",
    # orbits
    "compute_orbits",
    "dart_orbits",
    "vertex_orbits",
    "edge_orbits",
    "face_orbits",
    # classification
    "classify_symmetry",
    "UnknownSymmetryError",
    # point groups
    "PointGroup",
    "cyclic",
    "dihedral",
    "tetrahedral",
    "octahedral",
    "icosahedral",
]
