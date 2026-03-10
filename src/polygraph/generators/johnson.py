"""Johnson solid generators.

Parametric generators for pyramid and bipyramid families, plus individual
stubs for the Johnson solids that are convex deltahedra but have no simple
parametric form.  Stubs raise ``NotImplementedError`` until implemented.

Each function will return ``DartMap.from_face_lists(faces, num_vertices)``.
"""

from __future__ import annotations


def pyramid(n):
    """Return the n-gonal pyramid.  V=n+1, E=2n, F=n+1."""
    raise NotImplementedError


def bipyramid(n):
    """Return the n-gonal bipyramid.  V=n+2, E=3n, F=2n."""
    raise NotImplementedError


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
