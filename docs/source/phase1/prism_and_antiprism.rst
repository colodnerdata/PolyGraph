Phase 1a: Prism and Antiprism Families
======================================

Module: ``polygraph.generators.prisms``

This module provides parametric constructors for two regular polyhedral
families:

- ``prism(n)`` with counts ``(V, E, F) = (2n, 3n, n+2)``
- ``antiprism(n)`` with counts ``(V, E, F) = (2n, 4n, 2n+2)``

Both require ``n >= 3`` and return a
``polygraph.structures.dart_map.DartMap`` built from oriented face cycles.

.. automodule:: polygraph.generators.prisms
   :members:
   :undoc-members:
   :show-inheritance:
