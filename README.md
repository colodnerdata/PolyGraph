# PolyGraph

PolyGraph is a research-oriented Python library for studying polyhedral graphs using combinatorial maps (dart maps). It is currently narrowly focused on convex simply connected surface polyhedra.

It provides tools for topology, symmetry analysis, planar embeddings, and geometric realization of polyhedra.

The goal of PolyGraph is to separate combinatorial structure, graph algorithms, and geometric realization, enabling experimentation with planar graph drawing algorithms and polyhedral construction methods.

PolyGraph is designed as a computational playground for research in:

- polyhedral topology
- planar graph drawing
- combinatorial maps
- symmetry and automorphism groups
- geometric realization of polyhedra

---

Motivation

Many graph and geometry libraries mix topology, layout algorithms, and geometry into a single representation.

PolyGraph instead follows a layered architecture:

combinatorial topology → algorithms → geometric realization → visualization

This separation makes it easier to implement and compare different algorithms for:

- planar embeddings
- convex grid drawings
- symmetry detection
- polyhedral realizations

---

Core Concepts

Dart Maps  
A dart map (combinatorial map) represents the topology of a surface using directed half-edges.

Planar Embedding  
Planar graph drawing algorithms compute vertex coordinates while preserving the embedding.

Polyhedral Realization  
Geometric solvers recover vertex positions in 3D that satisfy face and symmetry constraints.

Symmetry Analysis  
Automorphism detection identifies vertex, edge, and face orbits.

---

Architecture

src/polygraph/

structures/
    Core combinatorial map representation

generators/
    Polyhedron generators

algorithms/
    Graph algorithms and symmetry analysis

geometry/
    Planar layout and 3D realization

visualization/
    Rendering utilities

export/
    Export to common mesh and graph formats

---

Example

```python
from polygraph.generators.platonic import cube

mesh = cube()

print(mesh.num_vertices())
print(mesh.num_faces())
```

---

Algorithms and Literature

PolyGraph aims to provide implementations of classical and modern planar graph drawing algorithms.

Canonical Ordering Methods
Incremental embeddings of 3-connected planar graphs based on canonical orderings.

Chrobak–Kant Convex Grid Drawing
Chrobak, M., & Kant, G. (1997). Convex grid drawings of 3-connected planar graphs. International Journal of Computational Geometry & Applications.

Disk-Link Convex Grid Drawing
Bekos, M. A., Gronemann, M., Montecchiani, F., & Symvonis, A.
Convex Grid Drawings of Planar Graphs with Constant Edge-Vertex Resolution.
Theoretical Computer Science 982:114290, 2024.
DOI: 10.1016/j.tcs.2023.114290

These algorithms provide a foundation for exploring:

- convex planar embeddings
- integer grid drawings
- symmetry-aware layouts
- polyhedral graph visualization

---

Research Roadmap

Near-term goals:

- DartMap topology implementation
- symmetry detection via automorphism groups
- planar embeddings using canonical ordering
- convex grid drawing algorithms

Future exploration:

- disk-link convex grid drawings
- force-directed planar refinement
- polyhedral realization via optimization
- symmetry-aware embeddings
- interactive algorithm visualization

---

Status

PolyGraph is an early-stage research project.

Interfaces may change rapidly while the architecture stabilizes.

---

License

Apache License 2.0