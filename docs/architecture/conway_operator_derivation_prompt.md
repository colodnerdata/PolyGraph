# Prompt: Combinatorial algorithms for the nine primitive Conway operators

This document is a self-contained brief for a mathematics student or another
LLM. It asks for a **derivation and specification** (not yet an implementation)
of the nine primitive Conway operators on combinatorial maps, with explicit
rules for *which new vertices are created* and *which raw constructions must be
identified (merged)*, including the consequences for edges and faces.

---

## 1. The task in one paragraph

PolyGraph represents a polyhedron as a `DartMap`, with no coordinates. A
Conway operator is naturally described geometrically ("put a vertex at each
edge midpoint", "cut each corner"). Translated to darts, such a description
produces **raw vertices that must be glued**: the same geometric vertex is
reached from several darts, edges or faces. Example: ambo's edge-midpoint
vertex is reached from both `d` and `alpha(d)`. For each operator, derive a
purely combinatorial rule that says:

1. **Create**: which new vertices exist, and what old object each is keyed on.
2. **Identify / merge**: which raw keys name the same new vertex, and which new
   vertices are just old vertices carried over.
3. **Edges and faces**: the new faces as cyclic lists of merged vertex ids
   (edges follow from faces), correctly oriented, so that each directed edge
   occurs exactly once and `DartMap.from_face_lists` accepts the result.
4. **Provenance**: for every new vertex, dart, edge and face, which old
   object(s) it came from (see section 5).

The algorithms must work for **every convex genus-0 polyhedron**, not only the
regular ones (section 7 explains why regular solids are a weak test set).

Operators in scope (`V, E, F` are the old counts; darts = `2E`):

| Op | Conway | Effect | V′ | E′ | F′ |
|---|---|---|---|---|---|
| dual | `d` | swap vertices and faces | F | E | V |
| ambo | `a` | vertices at edge midpoints | E | 2E | V+F |
| kis | `k` | apex on every face | V+F | 3E | 2E |
| truncate | `t` | cut every vertex | 2E | 3E | V+F |
| expand | `e` (= `aa`) | separate faces, fill gaps with quads | 2E | 4E | V+E+F |
| zip | `z` (= `dk`) | dual of kis | 2E | 3E | V+F |
| chamfer | `c` | every edge becomes a hexagon | V+2E | 4E | F+E |
| snub | `s` | chiral, gaps filled with triangles | 2E | 5E | 2+3E |
| whirl | `w` | chiral chamfer | V+4E | 7E | F+2E |

(Counts are from `docs/IMPLEMENTATION_ROADMAP.md`, Phase 10. They satisfy
V′−E′+F′ = 2 given V−E+F = 2. Re-derive them; do not trust them blindly.)

`kis` here means kis on **all** faces. The optional `n` parameter in the
skeleton doc (kis only n-gons) is out of scope. `gyro` and the algebraic
composites (`join`, `needle`, `ortho`, `meta`, `bevel`) are out of scope except
as cross-checks.

---

## 2. Concepts to settle before any formula

Clarify these in your answer before deriving anything.

**Darts, sigma, alpha, phi.** (See `docs/architecture/dart_map.md`.)

- A dart is an oriented half-edge, an integer `0..2E-1`. `alpha` pairs the two
  darts of an edge.
- `sigma(d)` is the next dart around `tail(d)`. `phi(d) = sigma^{-1}(alpha(d))`
  is the next dart along the face boundary.
- Vertices = cycles of `sigma`. Faces = cycles of `phi`. Edges = pairs
  `{d, alpha(d)}`.
- From `from_face_lists`: for a face `[.., t, u, v, ..]`, the dart `d = (u, v)`
  has `tail(d) = u`, `head(d) = v`, `phi(d)` = the dart `(v, w)`, and
  `sigma(d) = alpha(phi^{-1}(d)) = (u, t)`. The face containing `d` is the face
  **to the left of `d`**.
- A dart is therefore a *corner*: the triple (`tail(d)`, the edge, `face(d)`).
  Think of it as half of a flag. All locally-defined new vertices should be
  keyed on darts.

**Orientation convention (must be fixed first).** Chirality needs a global
orientation. Assume:

> Faces are listed counterclockwise when viewed from outside the polyhedron.
> Then `sigma` rotates counterclockwise about a vertex, and `phi` traverses a
> face counterclockwise.

The repository does not currently enforce or record this. Part of your work is
a short **audit**: for a convex realization of each generator (platonic,
prism, antiprism, pyramid, bipyramid), check by signed volume whether the face
lists are outward-CCW. Report any generator that is not. `dual_of`
(`structures/dual.py`) keeps the darts and sets `sigma' = phi`; establish and
state whether it preserves the orientation convention.

**Mirror.** `mirror(M)` has the same darts and `alpha`, with `sigma` replaced by
`sigma^{-1}` (equivalently, reverse every face list).

**Which operators are chiral?** Of the nine, only `snub` and `whirl`. The
other seven must satisfy `op(mirror(M)) ≅ mirror(op(M))`.

**Right-handed chirality (the user's choice).** Use right-handed for both
`snub` and `whirl`, consistently, and implement only that handedness. Define it
*operationally* so it can be falsified:

> Right-handed = in each old face, the inner "rotated copy" of the face is
> obtained by a rotation about the outward normal following the right-hand rule
> (counterclockwise seen from outside). Combinatorially, the new vertex
> belonging to dart `d` is displaced from `tail(d)` **toward `head(d)`**, inside
> `face(d)`.

State this definition in your answer and flag it for confirmation, since
"left" and "right" are named inconsistently in the literature. Left-handed
versions are obtained as `mirror ∘ op_R ∘ mirror`, never implemented
separately.

**Convex genus 0 means 3-connected planar.** By Steinitz's theorem, the
combinatorial types of convex polyhedra are exactly the 3-connected planar
graphs. So inputs are simple graphs (no loops or multi-edges), every vertex has
degree ≥ 3, every face has ≥ 3 sides and is a simple cycle, and two faces share
at most one edge. Your constructions may rely on this, and **the outputs must
again be 3-connected planar**.

---

## 3. Recommended method (a "raw keys, then identify" pipeline)

For each operator, structure the derivation as:

1. **Raw keys.** For each dart `d`, list the raw new vertices it emits, each
   named by a tuple, e.g. `("mid", d)`, `("corner", d)`, `("apex", face(d))`.
   Decide for each key what old object owns it.
2. **Identification relation.** List exactly which raw keys coincide, as an
   equivalence relation on keys (ambo: `("mid", d) ~ ("mid", alpha(d))`). State
   each class size and justify it from local geometry. Also list old vertices
   that are carried over (kis, chamfer, whirl).
3. **Canonical vertex keys.** Prefer choosing canonical representatives so no
   union-find is needed at run time (e.g. ambo key = `min(d, alpha(d))`, i.e.
   `edge_of_dart(d)`). Never rely on `vertex_orbits()` list order; label
   old vertices and faces by `min` dart of their orbit.
4. **Faces.** Give each new face as a cyclic list of canonical keys, one rule per
   face *type* (e.g. "face from old vertex v", "face from old edge e").
   Orientation must be outward-CCW.
5. **Edges.** Derive them from the faces and state each new edge's old origin
   (shrunk old edge, new spoke, new rim, ...). Check each undirected edge
   appears in exactly two faces, once per direction.
6. **Labelling.** Give a deterministic rule that turns canonical keys into
   integers `0..V′-1`, then call `from_face_lists`.
7. **Proof obligations.** Prove: identification classes have the claimed sizes;
   the face lists form a closed oriented 2-manifold; V′, E′, F′ match the
   table; the result is again simple and 3-connected (or say precisely when it
   would not be).

---

## 4. Deliverable: one "spec card" per operator

Use exactly these headings for each of the nine operators:

1. **Informal definition**: what it does to V, E, F, with face-size and
   vertex-degree effects.
2. **Raw keys and owners**: the table from step 1.
3. **Identifications**: the equivalence classes, and the carried-over old
   vertices.
4. **New faces**: one cyclic rule per face type, in terms of dart operations
   (`sigma`, `alpha`, `phi`, `tail`, `head`, `face`).
5. **New edges**: origin of each.
6. **Counts**: V′, E′, F′ plus the predicted *multiset of vertex degrees* and
   *multiset of face sizes* as a function of the old degree and face-size
   multisets.
7. **Chirality**: "achiral" or the right-handed rule above.
8. **Provenance maps**: see section 5.
9. **Worked example by hand** on a non-regular polyhedron (below), with the
   dart table written out.
10. **Python sketch**: `def op(dm: DartMap) -> ConwayResult`, with NumPy-style
    docstrings.

Hand-worked examples (get the dart tables from the repo generators; neither is
regular): `pyramid(4)` (one degree-4 vertex, one square, four triangles) and
`prism(3)`. Add one polyhedron from the random corpus in section 7 where
vertex degrees and face sizes vary.

Finish with an **interaction summary**: which operators share a key scheme
(e.g. `truncate`, `expand`, `snub` and `zip` all have `2E` vertices; explain
how their gluings differ), and a table of the cross-check identities in
section 6.

---

## 5. Provenance (required)

Each operator returns the new `DartMap` plus maps back to the old one.
Suggested shape (adjust if you find something better):

```python
@dataclass(frozen=True)
class Origin:
    kind: Literal["vertex", "edge", "face", "dart"]
    old_id: int              # canonical id: min dart of orbit, or the dart
    role: str = ""           # e.g. "mid", "apex", "corner", "spoke"

@dataclass(frozen=True)
class ConwayResult:
    dm: DartMap
    vertex_origin: list[Origin]          # new vertex -> old object
    dart_origin: list[Origin | None]     # new dart -> old object, None if new
    face_origin: list[Origin]            # new face -> old object
    merged_from: list[tuple]             # new vertex -> raw keys identified into it
```

Requirements: every new vertex, dart and face has exactly one origin; the
number of new objects per old object is a stated constant per role (e.g. ambo:
exactly 2 raw keys per vertex, 1 vertex per edge); and all maps are
**equivariant under relabelling** of the input darts. Provenance is what lets
later code (symmetry/orbits, iterated operators) reuse the result.

---

## 6. Checks the derivation must pass

Every claim must survive these tests, run on the corpus in section 7.

**Structural invariants of the output:** valid `DartMap`; V′, E′, F′ as
tabulated; connected; Euler characteristic 2; simple graph; min degree ≥ 3;
all faces ≥ 3 sides; 3-connected (e.g. `networkx.node_connectivity ≥ 3`);
predicted degree and face-size multisets.

**Label-invariance (equivariance).** For a random permutation `pi` of the
darts, `op(relabel(M, pi))` must be orientation-preserving isomorphic to
`op(M)`, and the provenance maps must transform accordingly. This catches
dart-index dependence **even when the polyhedron is symmetric**.

**Algebraic cross-checks, up to orientation-preserving isomorphism.** Verify
each one yourself before relying on it:

- `d(d(M)) = M`
- `a(d(M)) ≅ a(M)` (ambo is self-dual)
- `e(M) ≅ a(a(M))` (expand is primitive here; use this as a cross-check)
- `t(M) ≅ d(k(d(M)))`
- `z(M) ≅ d(k(M)) ≅ t(d(M))`
- `d(a(M))` = join, `d(a(a(M)))` = ortho (roadmap definitions)

**Chirality (`s`, `w` only).**

1. Orientation-preserving isomorphism class: `s(M)` and `mirror(s(M))` are
   isomorphic as *unoriented graphs* but, for suitable `M` (e.g. the cube), are
   **not** orientation-preserving isomorphic.
2. Independent left-handed rule: write a throwaway test-only `s_L_ref` using
   the opposite displacement (vertex of `d` toward `tail(d)`), coded
   separately from `s_R`. Require
   `s_L_ref(M) ≅ mirror(s_R(mirror(M)))` for every `M` in the corpus.
   (Defining `s_L` as `mirror ∘ s_R ∘ mirror` and testing against that same
   composite would be tautological, since `mirror` is an involution.)
3. **Geometric handedness test.** On a convex realization with outward-CCW
   faces, place each new vertex by the rule in section 2 (for `snub`, vertex of
   `d` = `(1-t)·tail + t·head` pulled slightly toward the face centre, small `t`).
   Every output face must still be outward-CCW (positive signed volume
   against the centroid). A wrong gluing or wrong handedness produces flipped
   triangles.

---

## 7. Test corpus: why regular solids are not enough

The five Platonic solids have large symmetry groups. A wrong dart pairing
(e.g. one that is off by a symmetry, or twisted the wrong way around a vertex)
often still yields a *valid map isomorphic to the right answer*, and
count-based checks cannot see a mixed-up gluing. Also all regular solids have
a single vertex degree and a single face size, which hides every bug that
confuses "degree of v" with "size of f" or a dart's two endpoints.

Require, in addition, a seeded generator of **irregular convex genus-0
polyhedra**:

1. **Simplicial:** the convex hull of `n` seeded random points on the sphere
   (`scipy.spatial.ConvexHull`). All faces are triangles, vertex degrees vary,
   and the automorphism group is trivial with overwhelming probability.
2. **Simple:** the dual of (1): all vertices degree 3, face sizes vary.
3. **Mixed:** start from (1) and delete random edges (merging two faces) as long
   as the graph stays 3-connected and min degree stays ≥ 3. This gives both
   varying degrees and varying face sizes.
4. **Structured irregular:** `pyramid(n)` for large `n` (one huge-degree
   vertex), `prism(n)`, `antiprism(n)`, `bipyramid(n)`, and operator results
   fed back in (e.g. `kis` on a single face region can be checked by hand).
5. Assert non-triviality of the corpus: record `|Aut|` (via the symmetry
   module / pynauty) and require that a good fraction has `|Aut| = 1`.

Sizes: small (n = 5–8, hand-checkable) and moderate (n ≈ 20–40, property
tests). Keep tests deterministic via seeds.

The corpus needs convexity only for the geometric handedness test; the
combinatorial tests need only 3-connectivity (Steinitz).

---

## 8. Constraints from the repository

- Python 3.11+, NumPy-convention docstrings, Ruff (79 columns, per `pyproject.toml`), no `print`.
- Build output only via `DartMap.from_face_lists(faces, num_vertices)`; it
  rejects faces with fewer than 3 vertices, repeated directed edges and
  missing opposite edges, so a wrong identification fails loudly.
- Target files: `src/polygraph/generators/conway.py` (currently empty) and
  `tests/generators/test_conway.py`. Tests first (TDD): write the invariant
  tests above, watch them fail, then implement.
- Do not add dependencies without a reason. NetworkX is an optional extra
  (`planar`), and SciPy is a core dependency.

## 9. Assumptions to confirm with the project owner

1. Right-handed means counterclockwise-seen-from-outside, with the new vertex of `d`
   displaced toward `head(d)` (section 2).
2. Outward-CCW is the global face orientation convention (audit may show some
   generators violate it).
3. `ortho = daa` and `expand = aa` are used only as cross-checks.
4. Provenance maps (section 5) are required for all nine operators.
