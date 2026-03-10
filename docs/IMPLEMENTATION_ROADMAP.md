# PolyGraph Implementation Roadmap

## Context

PolyGraph is a research library for polyhedral topology built on combinatorial maps. The **Structures layer is complete** (permutation, dart_map, traversal — ~1060 lines). Everything else — generators, algorithms, geometry, visualization, export, interop — exists as empty stubs. This plan sequences the implementation of all remaining layers, respecting dependency order and minimizing mathematical machinery at each stage.

**Guiding principles for a research package:**
- Implement the simplest correct version first; optimize later
- Prefer explicit loops over clever abstractions — readability matters more than DRY
- Each phase should produce testable, usable artifacts (not just plumbing)
- Numpy for array math; avoid custom linear algebra wrappers

---

## Phase 1: Generators — Platonic Solids & Parametric Families

**Why first:** Every subsequent layer needs concrete DartMaps to test against. Currently the only test fixture is a hand-built tetrahedron in test_traversal.py.

### 1a. `generators/platonic.py`
Each function returns `DartMap.from_face_lists(faces, n_vertices)`.

| Function | V | E | F | Face sizes | Math needed |
|---|---|---|---|---|---|
| `tetrahedron()` | 4 | 6 | 4 | all 3 | Enumerate 4 triangular faces |
| `cube()` | 8 | 12 | 6 | all 4 | Enumerate 6 quad faces with consistent orientation |
| `octahedron()` | 6 | 12 | 8 | all 3 | Dual of cube: 8 triangular faces |
| `dodecahedron()` | 20 | 30 | 12 | all 5 | 12 pentagonal faces (use known vertex adjacency) |
| `icosahedron()` | 12 | 30 | 20 | all 3 | 20 triangular faces |

**Minimum math:** Just hardcoded oriented face lists. Verify via Euler characteristic (V - E + F = 2) and genus (0).

### 1b. `generators/prisms.py`
| Function | V | E | F | Notes |
|---|---|---|---|---|
| `prism(n)` | 2n | 3n | n+2 | Two n-gon caps + n quads |
| `antiprism(n)` | 2n | 4n | 2n+2 | Two n-gon caps + 2n triangles |

**Math:** Parametric face generation from cyclic index arithmetic (mod n).

### 1c. `generators/johnson.py` (subset)
| Function | Notes |
|---|---|
| `pyramid(n)` | n-gon base + n triangles, V=n+1 |
| `bipyramid(n)` | 2n triangles, V=n+2 |

Defer `cupola(n)` and `rotunda()` — they need more complex face winding.

### Testing
- Euler characteristic = 2 for every generated solid
- Genus = 0
- Expected V, E, F counts
- Alpha is involution, sigma covers all darts
- Round-trip: `face_orbits()` recovers input faces (up to cyclic rotation)

### Files
- `src/polygraph/generators/platonic.py`
- `src/polygraph/generators/prisms.py`
- `src/polygraph/generators/johnson.py`
- `tests/generators/test_platonic.py`
- `tests/generators/test_prisms.py`

---

## Phase 2: Symmetry Detection

**Why now:** Symmetry detection depends only on the DartMap structure and the generators from Phase 1. Establishing the automorphism pipeline early enables symmetry-constrained optimization in later phases (geometry, layout refinement) and provides a clean validation target using well-understood Platonic solid symmetry groups.

### 2a. `interop/bliss_adapter.py`
BLISS (and nauty) compute automorphism generators via canonical labeling — use them directly rather than brute-forcing:
- **Option A (preferred):** Use `pynauty` or `python-bliss` bindings. Build a graph from the DartMap in a way that uniquely encodes the permutations `sigma` and `alpha`:
  - One vertex per dart.
  - Encode `sigma` and `alpha` as distinct edge types so that any graph automorphism necessarily commutes with both permutations.
  - If the backend supports directed and/or colored edges, use directed, edge-colored arcs: e.g., a directed edge of color `SIGMA` from `d` to `sigma(d)` and a directed edge of color `ALPHA` from `d` to `alpha(d)`.
  - If the backend only supports undirected, uncolored graphs, use a gadget/vertex-coloring scheme that breaks the symmetry between `sigma` vs `sigma^{-1}` and between `sigma` vs `alpha` (for example, introduce intermediate “sigma-edge” and “alpha-edge” vertices with distinct colors, or split each dart into in/out variants per permutation). The crucial requirement is that the encoding does **not** admit extra symmetries beyond DartMap automorphisms.
  - Call the library on this encoded graph to obtain automorphism generators and map them back to permutations on darts.
- **Option B (fallback only):** Pure-Python backtracking with pruning for environments where BLISS/nauty are unavailable (viable for polyhedra ≤ 120 darts), using the **same DartMap→graph encoding** as in Option A so that the computed automorphism group matches the BLISS/nauty semantics.

Target Option A from the start; Option B is a compatibility fallback, not the primary path. The well-understood Platonic solid examples (tetrahedron, cube, icosahedron) serve as validation against known group orders.

### 2b. `algorithms/symmetry/automorphisms.py`
- `compute_automorphism_generators(dm) -> list[Permutation]`: Return generating set for Aut(dm)
- `automorphism_group_order(generators, num_darts) -> int`: Compute |Aut(dm)| via Schreier-Sims or orbit-counting
- `is_orientation_preserving(pi, dm) -> bool`: Check if automorphism preserves face orientation

### 2c. `algorithms/symmetry/orbits.py`
- `compute_orbits(generators, n) -> list[list[int]]`: General orbit computation under a permutation group
- `dart_orbits(generators, dm)`, `vertex_orbits(generators, dm)`, `edge_orbits(generators, dm)`, `face_orbits(generators, dm)`: Derived orbit partitions

**Math:** Union-find on elements under generator action. For each generator g and element x, union x with g(x).

### Testing
- Tetrahedron: |Aut| = 24, 1 vertex orbit, 1 edge orbit, 1 face orbit
- Cube: |Aut| = 48, 1 vertex orbit, 1 edge orbit, 1 face orbit
- Prism(5): |Aut| = 20, 2 vertex orbits, 2 edge orbits, 2 face orbits

### Files
- `src/polygraph/interop/bliss_adapter.py`
- `src/polygraph/algorithms/symmetry/automorphisms.py`
- `src/polygraph/algorithms/symmetry/orbits.py`
- `tests/algorithms/test_symmetry.py`

---

## Phase 3: Symmetry Classification & Point Groups

### 3a. `algorithms/symmetry/point_groups.py`
- Define point group templates as named tuples: `(name, order, generators_pattern)`
- Families: C_n, D_n, T, O, I (and their variants with reflections)
- `cyclic(n)`, `dihedral(n)`, `tetrahedral()`, `octahedral()`, `icosahedral()`

### 3b. `algorithms/symmetry/classify.py`
- `classify_symmetry(generators, dm) -> str`: Match computed automorphism group against known point group templates
- Uses group order + orbit structure to narrow candidates

**Math:** Group order determines the family. |Aut|=24 on a solid with 1 vertex orbit → tetrahedral. |Aut|=48 → octahedral. |Aut|=120 → icosahedral. |Aut|=2n with 2 vertex orbits → dihedral D_n.

### Files
- `src/polygraph/algorithms/symmetry/point_groups.py`
- `src/polygraph/algorithms/symmetry/classify.py`
- `tests/algorithms/test_symmetry_classify.py`

---

## Phase 4: Structure Completions — Dual & Validation

### 4a. `structures/dual.py`
Dual map construction: swap vertex and face roles.

```
dual_sigma = phi        (face walk becomes vertex rotation)
dual_alpha = alpha      (edge involution unchanged)
dual_phi   = sigma      (vertex rotation becomes face walk)
```

- `dual_of(dm) -> DartMap`: Construct dual by reinterpreting permutations
- The dart set is the same; only sigma/alpha roles change

**Math:** For a combinatorial map (sigma, alpha), the dual is (phi, alpha) where phi = sigma⁻¹ ∘ alpha. Since our DartMap stores sigma and alpha, the dual's sigma is the original phi, and the dual's alpha stays the same.

**Validation:** `dual_of(cube())` should have the topology of an octahedron (V=6, E=12, F=8). `dual_of(dual_of(dm))` should recover the original.

### 4b. `structures/validation.py`
Extract and expose the invariant checks already embedded in `DartMap.__post_init__`:

- `check_alpha_involution(alpha)` — alpha[alpha[d]] == d, no fixed points
- `check_sigma_permutation(sigma)` — valid permutation on 0..n-1
- `check_euler_characteristic(dm, expected)` — V - E + F == expected
- `is_3_connected(dm)` — needed later for planar drawing (Chrobak-Kant requires 3-connectivity). Implement via checking that removing any 2 vertices leaves the graph connected. For small polyhedra this is tractable with BFS.

### Files
- `src/polygraph/structures/dual.py`
- `src/polygraph/structures/validation.py`
- `tests/structures/test_dual.py`
- `tests/structures/test_validation.py`

---

## Phase 5: Triangulation

**Why now:** Chrobak-Kant convex grid drawing (Phase 7: Canonical Ordering & Convex Grid Drawing) requires a **triangulated** 3-connected planar graph. Most Platonic solids are not triangulated (cube has quad faces, dodecahedron has pentagonal faces).

### 5a. `algorithms/triangulation/augment.py`
- `triangulate(dm) -> DartMap`: Fan-triangulate every non-triangular face by adding edges from one vertex to all non-adjacent vertices on the face boundary
- Track which edges/darts are "dummy" (added by triangulation) so they can be removed later — return a `TriangulationResult(dm, dummy_edges)` namedtuple

**Math:** For a face with k vertices (k > 3), pick one vertex v and add edges to the k-3 non-adjacent vertices. This adds k-3 edges and splits the face into k-2 triangles. New darts need correct sigma/alpha wiring.

**Implementation approach:** Rebuild face lists with subdivided faces, then call `DartMap.from_face_lists()`. This is simpler than in-place dart surgery and leverages the existing validated constructor.

### 5b. `algorithms/triangulation/validation.py`
- `is_triangulated(dm) -> bool`: Check that every face orbit has exactly 3 darts

### Testing
- `triangulate(cube())` should have F=12 (6 quads → 12 triangles), all faces triangular
- `triangulate(tetrahedron())` should be a no-op (already triangulated)
- Euler characteristic preserved
- `is_triangulated()` returns True after triangulation

### Files
- `src/polygraph/algorithms/triangulation/augment.py`
- `src/polygraph/algorithms/triangulation/validation.py`
- `tests/algorithms/test_triangulation.py`

---

## Phase 6: Planar Embedding Infrastructure

### 6a. `algorithms/planar/embedding.py`
- `PlanarEmbeddingView(dm)`: Thin wrapper providing neighbor-order and face-boundary queries in terms the planar algorithms expect. Mostly delegates to traversal functions.
- `ordered_neighbors(v)` → vertices adjacent to v in cyclic (sigma) order
- `face_boundary_vertices(f)` → boundary vertices of face f
- `degree(v)` → number of darts at vertex v

### 6b. `algorithms/planar/outer_face.py`
- `choose_outer_face(dm) -> int`: Select the face to serve as the outer (unbounded) face for planar drawing. Heuristic: choose the face with the most boundary vertices.
- `outer_face_anchors(dm, outer_face) -> tuple[int, int, int]`: Pick three vertices on the outer face boundary to anchor the convex drawing (needed by canonical ordering).

### Testing
- Verify cyclic neighbor order matches sigma orbit order
- Outer face of cube should be a 4-cycle
- Anchors are distinct vertices on the outer face boundary

### Files
- `src/polygraph/algorithms/planar/embedding.py`
- `src/polygraph/algorithms/planar/outer_face.py`
- `tests/algorithms/test_planar_embedding.py`

---

## Phase 7: Canonical Ordering & Convex Grid Drawing

This is the core near-term research target from the README.

### 7a. `algorithms/planar/canonical_order.py`
Canonical ordering of a **triangulated** 3-connected planar graph (de Fraysseix, Pach, Pollack / Kant).

**Input:** Triangulated DartMap + outer face with 3 anchor vertices (v1, v2, vn)
**Output:** Vertex ordering [v1, v2, ..., vn] such that:
  - v1, v2 are on the outer face and adjacent
  - vn is on the outer face
  - For k ≥ 3, vk has ≥ 2 neighbors in {v1, ..., v_{k-1}} and they form a contiguous interval on the outer boundary of G_{k-1}

**Algorithm:** Process vertices in reverse order. Maintain the outer boundary. At each step remove a vertex whose neighbors on the current outer boundary form a contiguous interval of length ≥ 2.

**Math:** Graph theory — contiguous neighbor intervals on a boundary path. Implementation uses the DartMap's face/vertex traversal to identify boundary structure.

### 7b. `algorithms/planar/contour_state.py`
- Track the evolving outer contour during incremental vertex placement
- Doubly-linked list of boundary vertices with efficient insert/delete
- `ContourState`: init from base edge (v1, v2), then `add_vertex(vk, left, right)` to update boundary

### 7c. `algorithms/planar/shift_structure.py`
- Shift accumulator for Chrobak-Kant x-coordinate computation
- Each vertex gets a relative x-shift; final coordinates computed by prefix sum

### 7d. `geometry/planar/layout.py`
- `chrobak_kant_layout(dm, outer_face=None) -> dict[int, tuple[int, int]]`
  - Full pipeline: triangulate → canonical order → shift-based placement
  - Returns integer grid coordinates for each vertex
  - All faces convex, grid size at most (2n-4) × (n-2)

**Math:** The Chrobak-Kant algorithm:
1. Compute canonical ordering
2. Place v1 at (0, 0), v2 at (2(n-1)-4, 0)
3. For each vk in order: place on lines of slope +1 and -1 from its leftmost and rightmost existing neighbors
4. Shift existing vertices to make room (shift structure)

### Testing
- All faces of the resulting drawing should be convex (cross-product test on consecutive edges)
- Coordinates are integers
- No edge crossings (for small cases, check all pairs)
- Works on triangulated tetrahedron, cube, octahedron, icosahedron

### Files
- `src/polygraph/algorithms/planar/canonical_order.py`
- `src/polygraph/algorithms/planar/contour_state.py`
- `src/polygraph/algorithms/planar/shift_structure.py`
- `src/polygraph/geometry/planar/layout.py`
- `tests/algorithms/test_canonical_order.py`
- `tests/geometry/test_planar_layout.py`

---

## Phase 8: Visualization (Planar)

### 8a. `visualization/matplotlib_planar.py`
- `draw_planar_graph(dm, positions, ax=None)`: Draw edges as line segments, vertices as points
- `draw_faces(dm, positions, ax=None, colors=None)`: Fill faces with color
- Optional: label vertices, highlight orbits

**Math:** None beyond coordinate indexing. Uses matplotlib patches/lines.

### 8b. `visualization/svg_planar.py`
- `render_svg(dm, positions) -> str`: Pure SVG string output (no matplotlib dependency)
- Edges as `<line>`, vertices as `<circle>`, faces as `<polygon>`

### Testing
- Smoke test: generates valid SVG/matplotlib figure without error
- Visual spot-check on cube, dodecahedron

### Files
- `src/polygraph/visualization/matplotlib_planar.py`
- `src/polygraph/visualization/svg_planar.py`
- `tests/visualization/test_planar_viz.py`

---

## Phase 9: Conway Operators

### 9a. `generators/conway.py`
Each operator takes a DartMap and produces a new DartMap.

| Operator | Effect | New V, E, F |
|---|---|---|
| `dual(dm)` | Swap vertices/faces | F, E, V |
| `kis(dm)` | Raise pyramids on faces | V+F, 3E, 2E (if triangulated) |
| `truncate(dm)` | Cut vertices | ... |
| `ambo(dm)` | Midpoints of edges become vertices | E, 2E, V+F |
| `expand(dm)` | ambo(ambo) | ... |
| `snub(dm)` | Expand + triangulate gaps | ... |

**Implementation:** Build new face lists from the old DartMap's combinatorial structure, then call `from_face_lists()`. The dual can delegate to `structures/dual.py`.

**Math:** Each operator is defined by a recipe that maps old darts/vertices/edges/faces to new face lists. The key insight is that all information is available from traversal — no geometry needed.

### Files
- `src/polygraph/generators/conway.py`
- `tests/generators/test_conway.py`

---

## Phase 10: Polyhedral Realization (3D Geometry)

### 10a. `geometry/polyhedral/face_planes.py`
- `FacePlaneParams`: Parameterize each face by a normal direction (2 spherical angles) + offset
- `expand_face_planes(params) -> tuple[ndarray, ndarray]`: Convert to (normals, offsets)
- `normal_from_spherical(phi, psi) -> ndarray`: Unit normal from angles

### 10b. `geometry/polyhedral/vertex_recovery.py`
- `recover_vertices(normals, offsets, dm) -> ndarray`: For each vertex, find the intersection of its incident face planes (least-squares solve of 3+ plane equations)

**Math:** Each vertex is the intersection of ≥ 3 planes. Solve n_i · x = d_i as a linear system. For exactly 3 planes: direct 3×3 solve. For >3: least-squares via `numpy.linalg.lstsq`.

### 10c. `geometry/polyhedral/initialization.py`
- Initialize face plane parameters from known geometry (e.g., Platonic solid coordinates)
- Or: random initialization with constraints (all normals pointing outward)

### 10d. `geometry/polyhedral/optimizer.py`
- `realize(dm, symmetry=None) -> ndarray`: Find 3D vertex positions satisfying:
  - Faces are planar
  - Edges have roughly uniform length (edge_uniformity_energy)
  - Dihedral angles are positive / faces don't intersect (dihedral_margin_energy)
  - Optional: symmetry constraints (reduce parameter space)
- Uses `scipy.optimize.minimize` with gradient-free or L-BFGS

**Math:** Energy minimization. Parameters: face plane angles + offsets. Objective: sum of edge-length variance + dihedral penalty. Vertex positions are derived quantities (from vertex_recovery).

### Files
- `src/polygraph/geometry/polyhedral/face_planes.py`
- `src/polygraph/geometry/polyhedral/vertex_recovery.py`
- `src/polygraph/geometry/polyhedral/initialization.py`
- `src/polygraph/geometry/polyhedral/optimizer.py`
- `tests/geometry/test_polyhedral.py`

---

## Phase 11: Layout Refinement & Advanced Drawing

### 11a. `geometry/planar/refinement.py`
- `refine_planar_layout(positions, dm, generators=None) -> ndarray`: Post-process grid drawing with force-directed smoothing
- `symmetry_energy(positions, generators)`: Penalize asymmetry
- `angular_resolution_energy(positions, dm)`: Maximize minimum angle at vertices

### 11b. `geometry/planar/layout.py` additions
- `disk_link_layout(dm)`: Bekos et al. algorithm (constant edge-vertex resolution)

### Files
- `src/polygraph/geometry/planar/refinement.py`
- `src/polygraph/geometry/planar/objectives.py`
- `src/polygraph/geometry/planar/constraints.py`

---

## Phase 12: Export & Interop

### 12a. Export formats
- `export/obj.py`: Wavefront OBJ for 3D meshes (vertices + face indices)
- `export/mesh_json.py`: JSON with vertices, faces, edges for web viewers
- `export/planar_json.py`: JSON with 2D positions + graph structure
- `export/svg.py`: Standalone SVG export

### 12b. Interop
- `interop/networkx_adapter.py`: Convert DartMap ↔ NetworkX graph (lose embedding info, keep adjacency)
- `structures/io.py`: Serialize/deserialize DartMap to JSON (face lists + vertex count)

### Files
- `src/polygraph/export/*.py`
- `src/polygraph/interop/networkx_adapter.py`
- `src/polygraph/structures/io.py`

---

## Phase 13: 3D Visualization

- `visualization/mesh_threejs.py`: Export three.js BufferGeometry JSON for web rendering
- Optional: `visualization/matplotlib_3d.py` for quick 3D plots via mpl_toolkits.mplot3d

---

## Verification Strategy

Each phase includes tests. The end-to-end pipeline test after Phase 8:

```python
from polygraph.generators.platonic import cube
from polygraph.algorithms.triangulation.augment import triangulate
from polygraph.geometry.planar.layout import chrobak_kant_layout
from polygraph.visualization.matplotlib_planar import draw_planar_graph

dm = cube()
tri_dm = triangulate(dm).dart_map
positions = chrobak_kant_layout(tri_dm)
draw_planar_graph(tri_dm, positions)
```

After Phase 10:
```python
from polygraph.geometry.polyhedral.optimizer import realize
vertices_3d = realize(cube())
```

Run at each phase:
```bash
python -m pytest tests/ -v
python -m ruff check src/ tests/
python -m mypy src/polygraph/
```
