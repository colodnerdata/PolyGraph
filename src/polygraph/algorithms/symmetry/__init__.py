"""Symmetry analysis for combinatorial maps.

Primary entry points
--------------------
:func:`~polygraph.algorithms.symmetry.automorphisms.compute_automorphism_generators`
    Compute a generating set for ``Aut(dm)`` via pynauty.
:func:`~polygraph.algorithms.symmetry.automorphisms.automorphism_group_order`
    Compute ``|Aut(dm)|`` from generators.
:func:`~polygraph.algorithms.symmetry.orbits.vertex_orbits`
    Symmetry orbits of vertices.
:func:`~polygraph.algorithms.symmetry.orbits.edge_orbits`
    Symmetry orbits of edges.
:func:`~polygraph.algorithms.symmetry.orbits.face_orbits`
    Symmetry orbits of faces.
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
