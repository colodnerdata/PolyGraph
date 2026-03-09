# PolyGraph Skeleton with Notional Methods

## Directory Structure

polygraph/\
├── pyproject.toml\
├── README.md\
├── OVERVIEW.md\
├── LICENSE\
│\
├── docs/\
│   └── source/\
│       ├── index.rst\
│       ├── api.rst\
│       ├── architecture.rst\
│       └── algorithms.rst\
│\
├── src/\
│   └── polygraph/\
│       ├── __init__.py\
│\
│       ├── structures/\
│       │   ├── __init__.py\
│       │   ├── dart_map.py\
│       │   ├── traversal.py\
│       │   ├── dual.py\
│       │   └── validation.py\
│       │\
│       ├── generators/\
│       │   ├── __init__.py\
│       │   ├── platonic.py\
│       │   ├── prisms.py\
│       │   ├── johnson.py\
│       │   └── conway.py\
│       │\
│       ├── algorithms/\
│       │   ├── __init__.py\
│       │   │\
│       │   ├── symmetry/\
│       │   │   ├── __init__.py\
│       │   │   ├── automorphisms.py\
│       │   │   ├── orbits.py\
│       │   │   ├── point_groups.py\
│       │   │   └── classify.py\
│       │   │\
│       │   └── planar/\
│       │       ├── __init__.py\
│       │       ├── embedding.py\
│       │       ├── outer_face.py\
│       │       └── canonical_order.py\
│       │\
│       ├── geometry/\
│       │   ├── **init**.py\
│       │   │\
│       │   ├── planar/\
│       │   │   ├── layout.py\
│       │   │   └── refinement.py\
│       │   │\
│       │   └── polyhedral/\
│       │       ├── face_planes.py\
│       │       ├── vertex_recovery.py\
│       │       └── optimizer.py\
│       │\
│       ├── visualization/\
│       │   ├── **init**.py\
│       │   ├── matplotlib_planar.py\
│       │   ├── mesh_threejs.py\
│       │   └── svg_planar.py\
│       │\
│       └── export/\
│           ├── **init**.py\
│           ├── obj.py\
│           └── mesh.py\
│\
└── tests/\
├── test_dart_map.py\
├── test_traversal.py\
├── test_dual.py\
├── test_generators.py\
├── test_symmetry.py\
├── test_planar_algorithms.py\
├── test_planar_layout.py\
└── test_polyhedral_realization.py\

---

# structures/

## dart_map.py

Core combinatorial map structure.

Notional methods:

* dart_map_from_face_lists(faces, num_vertices)
* edge_of_dart(d)
* vertex_orbits()
* face_orbits()
* euler_characteristic()
* genus()

Conceptual permutations:

sigma(d) → next dart around vertex
alpha(d) → opposite dart along edge
phi(d) → next dart around face

---

## traversal.py

Traversal utilities for combinatorial maps.

Notional methods:

* vertex_darts(dm, dart)
* face_darts(dm, dart)
* all_vertex_orbits(dm)
* all_face_orbits(dm)
* vertices_of_face(dm, face)
* faces_incident_to_vertex(dm, vertex)
* adjacent_face_pairs(dm)

---

## dual.py

Dual combinatorial map construction.

Notional methods:

* dual_of(dm)
* dual_sigma(dm)
* dual_phi(dm)
* dual_vertex_orbits(dm)
* dual_face_orbits(dm)

---

## validation.py

Topological invariant checks.

Notional methods:

* validate_dart_map(dm)
* check_alpha_involution(dm)
* check_sigma_permutation(dm)
* check_euler_characteristic(dm)
* is_3_connected(dm)

---

# generators/

## platonic.py

Platonic solid constructors.

Notional methods:

* tetrahedron()
* cube()
* octahedron()
* dodecahedron()
* icosahedron()

---

## prisms.py

Parametric prism families.

Notional methods:

* prism(n)
* antiprism(n)
* prism_faces(n)
* antiprism_faces(n)

---

## johnson.py

Johnson solid templates.

Notional methods:

* pyramid(n)
* cupola(n)
* rotunda()
* bipyramid(n)

---

## conway.py

Conway polyhedron operators.

Notional methods:

* dual(dm)
* ambo(dm)
* kis(dm, n=0)
* truncate(dm)
* expand(dm)
* snub(dm)
* bevel(dm)

---

# algorithms/symmetry/

## automorphisms.py

Notional methods:

* compute_automorphism_generators(dm)
* automorphism_group_order(generators, num_darts)
* is_orientation_preserving(pi, dm)

---

## orbits.py

Notional methods:

* compute_orbits(generators, elements)
* dart_orbits(generators, dm)
* vertex_orbits(generators, dm)
* edge_orbits(generators, dm)
* face_orbits(generators, dm)

---

## point_groups.py

Supported symmetry families:

Tetrahedral Symmetry (T, T_d, T_h)

* Symmetry of tetrahedron
* Group order: 24

Octahedral Symmetry (O, O_h)

* Symmetry of cube and octahedron
* Group order: 48

Icosahedral Symmetry (I, I_h)

* Symmetry of icosahedron and dodecahedron
* Group order: 120

Dihedral Symmetry (D_n, D_nd, D_nh)

* Symmetry of prisms and antiprisms
* Group order: 2n or 4n

Cyclic Symmetry (C_n, C_nv, C_nh, S_2n)

* Symmetry of pyramids
* Group order: n

Notional methods:

* cyclic(n, axis)
* dihedral(n, axis, perp)
* tetrahedral()
* octahedral()
* icosahedral()

---

## classify.py

Notional methods:

* classify_symmetry(generators, dm)
* concretize_symmetry(classification, face_orbit_reps, generators, dm)

---

# algorithms/planar/

## embedding.py

Notional methods:

* PlanarEmbeddingView(dm)
* ordered_neighbors(vertex)
* face_boundary_vertices(face)
* incident_faces(vertex)

---

## outer_face.py

Notional methods:

* choose_outer_face(dm, generators, embedding)
* score_outer_face(face)
* outer_face_anchors(face)

---

## canonical_order.py

Notional methods:

* canonical_order(dm, outer_face)
* build_shelling_trace(dm)
* canonical_steps(trace)

---

# geometry/planar/

## layout.py

Planar straight-line drawing algorithms.

Includes:

* Chrobak–Kant convex grid drawing
* Bekos et al. disk-link convex grid drawing

Notional methods:

* draw_planar(dm, outer_face=None)
* chrobak_kant_layout(dm)
* disk_link_layout(dm)
* initialize_contour_state(dm)

---

## refinement.py

Layout refinement and smoothing.

Notional methods:

* refine_planar_layout(positions, dm)
* symmetry_energy(positions, generators)
* angular_resolution_energy(positions, dm)

---

# geometry/polyhedral/

## face_planes.py

Plane parameter representation.

n · x = d

Notional methods:

* FacePlaneParams(theta)
* expand_face_planes(params)
* normal_from_spherical(phi, psi)
* make_rotation_matrix(omega)

---

## vertex_recovery.py

Recover vertices from plane intersections.

Notional methods:

* recover_vertices(normals, offsets, dm)
* recover_vertices_batch(normals, offsets)

---

## optimizer.py

Polyhedral realization optimizer.

Notional methods:

* realize(dm, symmetry=None)
* total_energy(params, dm)
* edge_uniformity_energy(positions, dm)
* dihedral_margin_energy(normals, dm)

---

# visualization/

## matplotlib.py

Notional methods:

* draw_planar_graph(dm, positions)
* draw_faces(dm, positions)
* draw_edges(dm, positions)

---

## threejs.py

Notional methods:

* to_threejs_geometry(mesh)
* export_threejs_scene(mesh)

---

# export/

## mesh.py

Mesh representation utilities.

Notional methods:

* to_mesh3d(dm, positions)
* triangulate(mesh)
* edge_list(mesh)

---

## obj.py

OBJ export utilities.

Notional methods:

* export_obj(mesh)
* write_obj_vertices(file, vertices)
* write_obj_faces(file, faces)

---

# Design Philosophy

PolyGraph separates combinatorial topology from geometric realization.

Layers:

1. Combinatorial structures
2. Graph algorithms
3. Planar drawing
4. Polyhedral realization
5. Visualization and export

This design supports experimentation with:

* planar graph drawing algorithms
* symmetry-aware geometry
* polyhedral realization methods
