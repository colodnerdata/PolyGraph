"""Triangulation algorithms for combinatorial dart maps."""

from __future__ import annotations

from polygraph.algorithms.triangulation.augment import (
    BarycentricResult,
    CellOrigin,
    CellType,
    barycentric_subdivision,
)
from polygraph.algorithms.triangulation.validation import (
    validate_triangulation,
)

__all__ = [
    "BarycentricResult",
    "CellOrigin",
    "CellType",
    "barycentric_subdivision",
    "validate_triangulation",
]
