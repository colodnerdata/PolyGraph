"""Public structures API."""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.dual import dual_of, is_isomorphism
from polygraph.structures.permutation import Permutation

__all__ = ["Permutation", "DartMap", "dual_of", "is_isomorphism"]
