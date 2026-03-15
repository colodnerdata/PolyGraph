"""Top-level PolyGraph public API."""
from __future__ import annotations

from polygraph.generators import (
    antiprism,
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    prism,
    tetrahedron,
)

__all__ = [
    "tetrahedron",
    "cube",
    "octahedron",
    "dodecahedron",
    "icosahedron",
    "prism",
    "antiprism",
]
