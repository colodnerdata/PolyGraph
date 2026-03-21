"""Public generators API."""
from __future__ import annotations

from polygraph.generators.johnson import bipyramid, dipyramid, pyramid
from polygraph.generators.notation import schlafli_symbol, vertex_configuration
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism

__all__ = [
    "tetrahedron",
    "cube",
    "octahedron",
    "dodecahedron",
    "icosahedron",
    "prism",
    "antiprism",
    "pyramid",
    "dipyramid",
    "bipyramid",
    "schlafli_symbol",
    "vertex_configuration",
]
