# Conway operator: dual (`d`)

Spec card for `polygraph.generators.conway.dual`, written to the format
requested in `conway_operator_derivation_prompt.md` so the other eight
primitive operators can follow the same template.

Source: `src/polygraph/generators/conway.py`.
Tests: `tests/generators/test_conway.py`, corpus in
`tests/generators/conway_corpus.py`.

---

## 1. Informal definition

The dual swaps the roles of vertices and faces. Put a new vertex in the middle
of every face, and join two new vertices by an edge when their faces share an
edge. Each old vertex is then surrounded by a new face.

| | V′ | E′ | F′ |
|---|---|---|---|
| dual | F | E | V |

A vertex of degree *k* becomes a *k*-gon, and a *k*-gon becomes a vertex of
degree *k*. Convexity (3-connectedness) is preserved: the dual of a convex
polyhedron is its polar polyhedron.

## 2. Raw keys and owners

The dual is the one primitive operator with **no raw vertices to reconcile**.
Every new vertex has exactly one key, and every new dart is an old dart.

| New object | Key | Owner | Count |
|---|---|---|---|
| vertex | `("face", f)` | old face `f`, canonical id `min` dart of its `phi`-orbit | F |
| dart | `d` | old dart `d` (same index) | 2E |
| face | `("vertex", v)` | old vertex `v`, canonical id `min` dart of its `sigma`-orbit | V |

## 3. Identifications

None. `ConwayResult.merged_from[i]` is the one-element tuple
`(("face", f),)` for every new vertex, and no old vertex is carried over
(unlike `kis`, `chamfer` and `whirl`).

## 4. The algorithm

Let `phi = sigma⁻¹ ∘ alpha`. The dual keeps the darts and `alpha` and uses
`phi` as the new vertex permutation:

```
(sigma*, alpha*) = (phi, alpha)
```

**New vertices.** The cycles of `sigma* = phi` are exactly the old faces. So
dart `d` now *starts* at the old face `face(d)` on its left, and the darts
leaving a new vertex appear in the counterclockwise boundary order of that
old face.

**New faces.** The dual's face permutation is

```
phi* = (sigma*)⁻¹ ∘ alpha = phi⁻¹ ∘ alpha = alpha ∘ sigma ∘ alpha
```

(using `phi⁻¹ = alpha⁻¹ ∘ sigma = alpha ∘ sigma`). This is `sigma`
conjugated by `alpha`, so the cycles of `phi*` are the sets `alpha(orbit)`
of old vertex orbits. The new face for the old vertex `v` consists of the
darts *pointing into* `v`: dart `d` has `head(d)` on its left.

```
        v = head(d)
           ^
           |  old dart d (south -> north)
   face(d) | face(alpha d)
   ------->|<-------        new dart d runs west -> east,
           |                crossing the old edge, with v on its left
        u = tail(d)
```

**New edges.** Old edge `{d, alpha(d)}` becomes the new edge with the same
two darts. It crosses the old edge, joining `face(d)` to `face(alpha(d))`.

**Pseudocode** (`O(n)` in the number of darts):

```
dual(dm):
    dm* = DartMap(sigma = [phi(d) for d in darts], alpha = dm.alpha)
    for each vertex orbit O of dm*:        # = an old face
        origin = ("face", min(O))
    for each face orbit O of dm*:          # = an old vertex
        origin = ("vertex", canonical id of the vertex containing alpha(O[0]))
    for each dart d:
        origin = ("dart", d)
```

Any one dart of a dual face identifies its old vertex, because all darts of a
dual face point into the same old vertex.

## 5. Orientation and chirality

Around an old face the darts leaving the new vertex appear in counterclockwise
order, so `sigma*` still turns counterclockwise and the dual keeps the
outward-counterclockwise convention. The dual is therefore **orientation
preserving** and **achiral**:

```
dual(mirror(M))  is  mirror(dual(M))   up to relabelling by alpha.
```

This is checked geometrically in the tests: a convex polyhedron is oriented
outward, dualised, and realised by polarity about an interior point, and every
dual face must again be counterclockwise from outside. A negative control
(the mirrored dual) must fail the same check.

## 6. Counts and invariants

- `(V′, E′, F′) = (F, E, V)`; `V′ − E′ + F′ = 2`.
- The multiset of dual vertex degrees equals the multiset of old face sizes,
  and the multiset of dual face sizes equals the multiset of old vertex
  degrees.

## 7. `dual(dual(M))` is not `M`

`dual(dual(M))` has `sigma = alpha ∘ sigma ∘ alpha`, so it equals `M`
**conjugated by `alpha`**, not `M` itself. Code that compares the two maps
dart by dart will fail; use

```python
is_isomorphism(dm, dual(dual(dm).dm).dm, dm.alpha)
```

## 8. Provenance maps

```python
result = dual(dm)
result.vertex_origin[i]  # Origin("face", f, "dual"), f = min dart of old face
result.face_origin[j]    # Origin("vertex", v, "dual"), v = min dart of old vertex
result.dart_origin[d]    # Origin("dart", d, "dual")
result.merged_from[i]    # (("face", f),)
```

Indices follow `result.dm.vertex_orbits()` and `result.dm.face_orbits()`.
Canonical ids are smallest darts, not list positions, so they do not depend on
the order orbits are listed.

Canonical ids do depend on dart *labels*. Under a relabelling `pi` of the
input darts, the **dart sets** named by the provenance move by `pi` (this is
tested), but the ids themselves do not.

## 9. Worked example: `pyramid(4)`

Faces (outward-CCW as listed): base `[0,1,2,3]` and sides `[1,0,4]`,
`[2,1,4]`, `[3,2,4]`, `[0,3,4]`. It has one degree-4 vertex (the apex 4), four
degree-3 vertices, one square and four triangles. The shape is self-dual, but
the dual's vertices and faces are *different objects* with different labels,
which is what the table below tracks.

| d | tail → head | face | `sigma` | `phi` | `alpha` | `sigma*` (=`phi`) | `phi*` |
|---|---|---|---|---|---|---|---|
| 0 | 0→1 | 0 | 13 | 1 | 4 | 1 | 6 |
| 1 | 1→2 | 0 | 4 | 2 | 7 | 2 | 9 |
| 2 | 2→3 | 0 | 7 | 3 | 10 | 3 | 12 |
| 3 | 3→0 | 0 | 10 | 0 | 13 | 0 | 15 |
| 4 | 1→0 | 1 | 8 | 5 | 0 | 5 | 3 |
| 5 | 0→4 | 1 | 0 | 6 | 15 | 6 | 14 |
| 6 | 4→1 | 1 | 15 | 4 | 8 | 4 | 7 |
| 7 | 2→1 | 2 | 11 | 8 | 1 | 8 | 0 |
| 8 | 1→4 | 2 | 1 | 9 | 6 | 9 | 5 |
| 9 | 4→2 | 2 | 6 | 7 | 11 | 7 | 10 |
| 10 | 3→2 | 3 | 14 | 11 | 2 | 11 | 1 |
| 11 | 2→4 | 3 | 2 | 12 | 9 | 12 | 8 |
| 12 | 4→3 | 3 | 9 | 10 | 14 | 10 | 13 |
| 13 | 0→3 | 4 | 5 | 14 | 3 | 14 | 2 |
| 14 | 3→4 | 4 | 3 | 15 | 12 | 15 | 11 |
| 15 | 4→0 | 4 | 12 | 13 | 5 | 13 | 4 |

- Old vertex orbits (cycles of `sigma`): `[0,13,5]`, `[1,4,8]`, `[2,7,11]`,
  `[3,10,14]`, `[6,15,12,9]` (the apex).
- Old face orbits (cycles of `phi`): `[0,1,2,3]` (base), `[4,5,6]`, `[7,8,9]`,
  `[10,11,12]`, `[13,14,15]`.
- **Dual vertices** (cycles of `sigma* = phi`) are exactly the old face orbits.
  Dual vertex 0 is the base, with origin `("face", 0)`, and has degree 4.
- **Dual faces** (cycles of `phi*`): `[0,6,7]`, `[1,9,10]`, `[2,12,13]`,
  `[3,15,4]`, `[5,14,11,8]`. The last one is the old apex 4: the darts
  pointing into 4 are 5, 8, 11, 14, and it is a square. The face `[3,15,4]`
  is the old vertex 0: darts 3 (3→0), 4 (1→0), 15 (4→0) point into it.
- Counts: `(5, 8, 5)` before and after, since a square pyramid is self-dual;
  the degree-4 apex became the square face and the square base became the
  degree-4 vertex.

## 10. Python

```python
from polygraph.generators.conway import dual
from polygraph.generators.johnson import pyramid

result = dual(pyramid(4))
result.dm                 # the dual DartMap
result.vertex_origin[0]   # Origin(kind='face', old_id=0, role='dual')
```

`dual` is built on `structures.dual.dual_of`, which does the dart-level
permutation swap; `conway.dual` adds the provenance.

## 11. Testing

Run on the five Platonic solids, prism/antiprism/pyramid/dipyramid families,
and 14 seeded irregular polyhedra (mixed vertex degrees up to 9 and face
sizes up to 8; more than half with trivial symmetry group).

| Check | Why it is not satisfied by symmetric inputs |
|---|---|
| counts, Euler, connected | basic invariants |
| degree and face-size multisets swap | catches vertex/face mix-ups |
| dual vertices are old face dart sets; dual faces are `alpha`(old vertex) | pins the exact gluing, not only isomorphism class |
| `dual(dual(M))` isomorphic to `M` via `alpha` | involution up to relabelling |
| commutes with `mirror` | achirality |
| relabel input darts by random `pi`, output is isomorphic via `pi` | dart-numbering independence, valid on symmetric solids |
| provenance dart sets move with `pi` | provenance is equivariant |
| geometric outward-CCW check, plus mirrored negative control | the only check that sees orientation |

The checks were mutation-tested: a dual built with the wrong chirality,
face provenance that forgets `alpha`, an off-by-one dart origin, and `max`
instead of `min` canonical ids each fail the suite.

## 12. Open items

- The generators (`platonic`, `prisms`, `johnson`) have not been audited for
  the outward-counterclockwise convention. The geometric test above builds its
  own oriented face lists, so it does not depend on them.
- `kis(dm, n)` and the composite operators (`join = d∘a`, `needle`, ...) can
  now build on `dual` and its provenance.
