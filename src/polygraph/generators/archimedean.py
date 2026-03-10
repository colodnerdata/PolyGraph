"""Archimedean solid generators.

All 13 Archimedean solids are vertex-transitive with regular polygon faces.
Several arise naturally from Conway operators on Platonic solids (truncate,
ambo, expand, snub — see generators/conway.py), so hardcoded implementations
are deferred until Phase 9 (Conway operators) is complete.

Each function will return ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations


def truncated_tetrahedron():
    """Return the truncated tetrahedron (3.6.6). V=12, E=18, F=8."""
    raise NotImplementedError


def cuboctahedron():
    """Return the cuboctahedron (3.4.3.4). V=12, E=24, F=14."""
    raise NotImplementedError


def truncated_cube():
    """Return the truncated cube (3.8.8). V=24, E=36, F=14."""
    raise NotImplementedError


def truncated_octahedron():
    """Return the truncated octahedron (4.6.6). V=24, E=36, F=14."""
    raise NotImplementedError


def rhombicuboctahedron():
    """Return the rhombicuboctahedron (3.4.4.4). V=24, E=48, F=26."""
    raise NotImplementedError


def truncated_cuboctahedron():
    """Return the truncated cuboctahedron (4.6.8). V=48, E=72, F=26."""
    raise NotImplementedError


def snub_cube():
    """Return the snub cube (3.3.3.3.4). V=24, E=60, F=38."""
    raise NotImplementedError


def icosidodecahedron():
    """Return the icosidodecahedron (3.5.3.5). V=30, E=60, F=32."""
    raise NotImplementedError


def truncated_dodecahedron():
    """Return the truncated dodecahedron (3.10.10). V=60, E=90, F=32."""
    raise NotImplementedError


def truncated_icosahedron():
    """Return the truncated icosahedron (5.6.6). V=60, E=90, F=32."""
    raise NotImplementedError


def rhombicosidodecahedron():
    """Return the rhombicosidodecahedron (3.4.5.4). V=60, E=120, F=62."""
    raise NotImplementedError


def truncated_icosidodecahedron():
    """Return the truncated icosidodecahedron (4.6.10). V=120, E=180, F=62."""
    raise NotImplementedError


def snub_dodecahedron():
    """Return the snub dodecahedron (3.3.3.3.5). V=60, E=150, F=92."""
    raise NotImplementedError
