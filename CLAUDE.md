# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install for development
python -m pip install -e .

# Run all tests
python -m pytest tests/ -v

# Run a single test file
python -m pytest tests/structures/test_traversal.py -v

# Run a single test by name
python -m pytest tests/structures/test_traversal.py::test_orbit_yields_full_cycle_without_repeating_start -v

# Lint
python -m ruff check src/ tests/

# Type check
python -m mypy src/polygraph/

# Build distribution
python -m build
```

## Architecture

PolyGraph is a research library for polyhedral topology built on **combinatorial maps** (dart maps). It enforces strict separation of concerns across four layers:

```
Visualization  (matplotlib, SVG, three.js)
     ↑
Geometry       (planar layout, 3D polyhedral realization)
     ↑
Algorithms     (symmetry, planar embedding, triangulation)
     ↑
Structures     (DartMap, Permutation, Traversal)  ← current focus
     ↑
Generators     (Platonic solids, prisms, Conway operators)
```

### Core Topology Model (`src/polygraph/structures/`)

The foundation. Everything else builds on these three types:

- **`permutation.py`** — Immutable finite permutation on integer indices. Supports `inverse()`, `compose()`, `orbit(start)`, `cycles()`. Used everywhere as the primitive operation type.

- **`dart_map.py`** — The central data structure. A combinatorial map where:
  - `sigma[d]` = next dart around the same vertex (cyclic vertex order)
  - `alpha[d]` = opposite dart on the same edge (fixed-point-free involution)
  - `phi(d) = sigma⁻¹(alpha(d))` = next dart around the same face (derived, not stored)

  Constructed via `DartMap.from_face_lists(faces, num_vertices)`. Supports orbit queries, Euler characteristic, and genus computation.

- **`traversal.py`** — Low-level iteration helpers using a **representative-dart convention** (each orbit is identified by its smallest dart index). Key functions: `orbit()`, `vertex_darts()`, `face_darts()`, `edge_darts()`, and enumerators `all_vertex_orbits()`, `all_face_orbits()`, `all_edge_orbits()`.

### Algorithms (`src/polygraph/algorithms/`)

- **symmetry/**: Automorphism group analysis, orbit classification, point group detection. Uses BLISS (via interop) for automorphism computation.
- **planar/**: Canonical ordering + Chrobak–Kant / Bekos et al. convex grid drawing for planar embeddings.
- **triangulation/**: Augmentation with dummy vertices/edges to produce triangulated maps.

### Linting Rules

Ruff is configured with rules `E, F, I, B, UP, N, D` and numpy docstring convention (line length 79). Docstring rules D203, D213, D100, D104 are ignored. D1xx rules for public members are **currently relaxed** (`per-file-ignores` for `structures/` modules) while the core layer is being built out — restore them when modules are complete.
