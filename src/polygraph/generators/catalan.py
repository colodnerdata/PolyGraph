"""Catalan solid generators.

The 13 Catalan solids are the duals of the 13 Archimedean solids. Hardcoded
implementations are deferred until Phase 2 dual construction (``dual_of``)
and Phase 10 Conway operators are complete; most can then be derived via
``dual_of(archimedean_solid())``.

Each function will return ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations


def triakis_tetrahedron():
    """Return the triakis tetrahedron (dual of truncated tetrahedron)."""
    raise NotImplementedError


def rhombic_dodecahedron():
    """Return the rhombic dodecahedron (dual of cuboctahedron)."""
    raise NotImplementedError


def triakis_octahedron():
    """Return the triakis octahedron (dual of truncated cube)."""
    raise NotImplementedError


def tetrakis_hexahedron():
    """Return the tetrakis hexahedron (dual of truncated octahedron)."""
    raise NotImplementedError


def deltoidal_icositetrahedron():
    """Return the deltoidal icositetrahedron (dual of rhombicuboctahedron)."""
    raise NotImplementedError


def disdyakis_dodecahedron():
    """Return the disdyakis dodecahedron (dual of truncated cuboctahedron)."""
    raise NotImplementedError


def pentagonal_icositetrahedron():
    """Return the pentagonal icositetrahedron (dual of snub cube)."""
    raise NotImplementedError


def rhombic_triacontahedron():
    """Return the rhombic triacontahedron (dual of icosidodecahedron)."""
    raise NotImplementedError


def triakis_icosahedron():
    """Return the triakis icosahedron (dual of truncated dodecahedron)."""
    raise NotImplementedError


def pentakis_dodecahedron():
    """Return the pentakis dodecahedron (dual of truncated icosahedron)."""
    raise NotImplementedError


def deltoidal_hexacontahedron():
    """Return the deltoidal hexacontahedron.

    Dual of the rhombicosidodecahedron.
    """
    raise NotImplementedError


def disdyakis_triacontahedron():
    """Return the disdyakis triacontahedron.

    Dual of the truncated icosidodecahedron.
    """
    raise NotImplementedError


def pentagonal_hexacontahedron():
    """Return the pentagonal hexacontahedron (dual of snub dodecahedron)."""
    raise NotImplementedError
