# Dart Map Architecture

## Purpose

`DartMap` is the core combinatorial data structure used by PolyGraph to represent embedded graphs and polyhedral topology without geometry.

It stores topology as permutations on **darts** (oriented half-edges), which makes vertex/edge/face traversal and topological invariant computation straightforward.

## Core Objects

- **Dart**: an oriented half-edge, identified by an integer index.
- **`sigma`**: permutation of darts; `sigma[d]` is the next dart around the same vertex.
- **`alpha`**: fixed-point-free involution; `alpha[d]` is the opposite dart on the same undirected edge.
- **`phi`**: derived face permutation in this implementation:
  - `phi(d) = sigma^{-1}(alpha(d))`
  - exposed as `DartMap.phi(d)` (computed, not stored)

## Structural Invariants

`DartMap.__post_init__` enforces:

- `sigma` and `alpha` have equal non-zero length.
- `sigma` is a valid permutation of `0..n-1`.
- `alpha` is a valid permutation and:
  - has no fixed points (`alpha[d] != d`)
  - is an involution (`alpha[alpha[d]] == d`)
- Number of darts is even (`n % 2 == 0`).

These conditions imply each undirected edge is represented by exactly two darts.

## Orbits and Topology

- **Vertex orbits** are cycles of `sigma` (`vertex_orbits()`).
- **Face orbits** are cycles of `phi` (`face_orbits()`).
- **Edges** are identified canonically by `min(d, alpha[d])` (`edge_of_dart()`).
- **Connected components** are discovered by graph traversal over `{sigma, sigma^{-1}, alpha}`.
- **Euler characteristic**:
  - `chi = V - E + F`
  - implemented by `euler_characteristic()`
- **Orientable genus**:
  - `g = (2c - chi) / 2`, where `c` is number of connected components
  - implemented by `genus()`

## Construction from Face Lists

`dart_map_from_face_lists(faces, num_vertices)` builds a valid `DartMap` from oriented face boundaries:

1. For each directed face edge `(u, v)`, allocate one dart.
2. Build face successor permutation `phi` from cyclic order within each face.
3. Build `alpha` by matching each `(u, v)` to `(v, u)`.
4. Recover `sigma` from:
   - `sigma = alpha ∘ phi^{-1}`
5. Return `DartMap(sigma, alpha)`, which re-validates all invariants.

Input assumptions:

- Faces have length at least 3.
- Vertex indices are in `0..num_vertices-1`.
- No repeated directed edge `(u, v)`.
- Every directed edge has opposite direction present, i.e. closed 2-manifold style input.

## Public API Surface

Class methods/properties:

- `num_darts`
- `num_edges`
- `phi(d)`
- `edge_of_dart(d)`
- `vertex_orbits()`
- `face_orbits()`
- `connected_components()`
- `euler_characteristic()`
- `genus()`
- `from_face_lists(faces, num_vertices)`

Module-level wrappers:

- `dart_map_from_face_lists(...)`
- `edge_of_dart(dm, d)`
- `vertex_orbits(dm)`
- `face_orbits(dm)`
- `euler_characteristic(dm)`
- `genus(dm)`

## Minimal Example

```python
from polygraph.structures.dart_map import DartMap

# Tetrahedron faces
faces = [
    [0, 1, 2],
    [0, 3, 1],
    [1, 3, 2],
    [0, 2, 3],
]

dm = DartMap.from_face_lists(faces, num_vertices=4)

assert dm.num_darts == 12
assert dm.num_edges == 6
assert len(dm.vertex_orbits()) == 4
assert len(dm.face_orbits()) == 4
assert dm.euler_characteristic() == 2
assert dm.genus() == 0
```

## Why This Representation

Permutation-based combinatorial maps separate connectivity from coordinates. That keeps this layer reusable across:

- planar embedding and traversal algorithms
- dual construction and symmetry analysis
- 3D realization and export pipelines

Geometry can be added later without changing the topological core.
