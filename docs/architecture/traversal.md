# Traversal Architecture

## Purpose

`polygraph.structures.traversal` provides low-level iteration helpers over
`DartMap` orbits and incidence relations.

The module sits in the `structures` layer because it only depends on
combinatorial permutations (`sigma`, `alpha`, `phi`) and not on geometry or
higher-level algorithms.

## Representative-Dart Model

Traversal APIs in this module treat vertices, faces, and edges as represented
by dart indices:

- vertex handle: any dart in a `sigma` orbit
- face handle: any dart in a `phi` orbit
- edge handle: any dart in an `alpha` pair

Handles are valid if they are in `0..dm.num_darts-1`.

## Core Traversal Semantics

- `orbit(start, perm)`:
  - yields a full cycle under `perm`, starting at `start`
  - stops before repeating `start`
- `vertex_darts(dm, d)`:
  - yields the `sigma` orbit of `d`
- `face_darts(dm, d)`:
  - yields the `phi` orbit of `d`
- `edge_darts(dm, d)`:
  - yields `d`, then `alpha[d]`

## Orbit Enumerators

- `all_vertex_orbits(dm)`:
  - yields one representative dart per vertex orbit
  - representatives are discovered by scanning darts in increasing order
- `all_face_orbits(dm)`:
  - same policy for face orbits
- `all_edge_orbits(dm)`:
  - canonical representative policy: yield `d` iff `d < alpha[d]`

These representatives are deterministic for a fixed labeling, but they are not
canonical under relabeling.

## Incidence/Adjacency Queries

- `neighbors(dm, vertex)`:
  - yields `alpha[d]` for each `d` in the source vertex orbit
  - preserves cyclic order and multiplicity (parallel edges repeat)
- `vertices_of_face(dm, face)`:
  - yields boundary vertex representatives in face-walk order
  - under this model, equivalent to `face_darts(dm, face)`
- `faces_incident_to_vertex(dm, vertex)`:
  - yields incident face representatives in rotational order
  - under this model, equivalent to `vertex_darts(dm, vertex)`
- `adjacent_face_pairs(dm)`:
  - yields one `(face_a, face_b)` pair per edge as `(e, alpha[e])`
  - where `e` is from `all_edge_orbits(dm)`

## Relationship to `DartMap`

`DartMap` defines:

- `sigma` for vertex rotation
- `alpha` for opposite darts on an edge
- `phi(d) = sigma^{-1}(alpha(d))` for face traversal

Traversal utilities provide a focused API for consuming these permutations
without duplicating topology logic in algorithm modules.

## Minimal Example

```python
from polygraph.structures.dart_map import DartMap
from polygraph.structures.traversal import (
    all_edge_orbits,
    all_face_orbits,
    all_vertex_orbits,
    neighbors,
)

faces = [
    [0, 1, 2],
    [0, 3, 1],
    [1, 3, 2],
    [0, 2, 3],
]
dm = DartMap.from_face_lists(faces, 4)

v_reps = list(all_vertex_orbits(dm))
f_reps = list(all_face_orbits(dm))
e_reps = list(all_edge_orbits(dm))

assert len(v_reps) == 4
assert len(f_reps) == 4
assert len(e_reps) == 6
assert len(list(neighbors(dm, v_reps[0]))) == 3
```
