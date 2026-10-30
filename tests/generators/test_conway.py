"""Tests for the Conway operators in :mod:`polygraph.generators.conway`.

Each operator is tested on regular solids *and* on an irregular corpus
(``conway_corpus``): symmetric inputs can hide wrong gluings, so every test
that can fail silently on a symmetric solid also runs on asymmetric ones.
"""

from __future__ import annotations

import random

import pytest

from polygraph.generators.conway import ConwayResult, Origin, dual
from polygraph.generators.johnson import dipyramid, pyramid
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism
from polygraph.structures.dart_map import DartMap
from polygraph.structures.dual import is_isomorphism
from tests.generators.conway_corpus import (
    irregular_corpus,
    irregular_polyhedron,
    mirror,
    orientation_preserving_aut_order,
    relabel,
)

REGULAR_AND_FAMILIES = [
    ("tetrahedron", tetrahedron),
    ("cube", cube),
    ("octahedron", octahedron),
    ("dodecahedron", dodecahedron),
    ("icosahedron", icosahedron),
    ("prism(3)", lambda: prism(3)),
    ("prism(5)", lambda: prism(5)),
    ("antiprism(4)", lambda: antiprism(4)),
    ("pyramid(4)", lambda: pyramid(4)),
    ("pyramid(7)", lambda: pyramid(7)),
    ("dipyramid(5)", lambda: dipyramid(5)),
]

ALL_MAPS = [(name, factory()) for name, factory in REGULAR_AND_FAMILIES] + (
    irregular_corpus()
)


@pytest.fixture(params=ALL_MAPS, ids=[name for name, _ in ALL_MAPS])
def any_dm(request) -> DartMap:
    return request.param[1]


@pytest.fixture(
    params=irregular_corpus(), ids=[name for name, _ in irregular_corpus()]
)
def irregular_dm(request) -> DartMap:
    return request.param[1]


def _counts(dm: DartMap) -> tuple[int, int, int]:
    return len(dm.vertex_orbits()), dm.num_edges, len(dm.face_orbits())


# ---------------------------------------------------------------------------
# Corpus sanity: the irregular inputs really are irregular.
# ---------------------------------------------------------------------------


def test_irregular_corpus_is_irregular() -> None:
    corpus = irregular_corpus()
    degree_sets = [
        {len(o) for o in dm.vertex_orbits()} for _, dm in corpus
    ]
    face_size_sets = [{len(o) for o in dm.face_orbits()} for _, dm in corpus]

    assert all(len(s) >= 2 for s in degree_sets)
    assert any(len(s) >= 2 for s in face_size_sets)
    asymmetric = [
        dm for _, dm in corpus if orientation_preserving_aut_order(dm) == 1
    ]
    assert len(asymmetric) >= len(corpus) // 2


def test_irregular_polyhedron_is_deterministic() -> None:
    a = irregular_polyhedron(3, 9, 2)
    b = irregular_polyhedron(3, 9, 2)
    assert (a.sigma, a.alpha) == (b.sigma, b.alpha)


# ---------------------------------------------------------------------------
# dual
# ---------------------------------------------------------------------------


def test_dual_returns_conway_result(any_dm) -> None:
    result = dual(any_dm)
    assert isinstance(result, ConwayResult)
    assert isinstance(result.dm, DartMap)


def test_dual_swaps_counts(any_dm) -> None:
    v, e, f = _counts(any_dm)
    assert _counts(dual(any_dm).dm) == (f, e, v)


def test_dual_is_a_sphere(any_dm) -> None:
    dm = dual(any_dm).dm
    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0
    assert len(dm.connected_components()) == 1


def test_dual_swaps_degree_and_face_size_multisets(any_dm) -> None:
    dm = dual(any_dm).dm
    assert sorted(len(o) for o in dm.vertex_orbits()) == sorted(
        len(o) for o in any_dm.face_orbits()
    )
    assert sorted(len(o) for o in dm.face_orbits()) == sorted(
        len(o) for o in any_dm.vertex_orbits()
    )


def test_dual_vertices_are_exactly_the_old_faces(any_dm) -> None:
    # sigma* = phi, so dual vertex orbits are the old face dart sets.
    result = dual(any_dm)
    assert {frozenset(o) for o in result.dm.vertex_orbits()} == {
        frozenset(o) for o in any_dm.face_orbits()
    }


def test_dual_faces_are_alpha_images_of_old_vertices(any_dm) -> None:
    # phi* = alpha sigma alpha, so the dual face of old vertex v is the set
    # of darts pointing *into* v.
    result = dual(any_dm)
    expected = {
        frozenset(any_dm.alpha[d] for d in orbit)
        for orbit in any_dm.vertex_orbits()
    }
    assert {frozenset(o) for o in result.dm.face_orbits()} == expected


def test_dual_darts_cross_the_same_edges(any_dm) -> None:
    result = dual(any_dm)
    assert result.dm.alpha == any_dm.alpha
    assert result.dm.num_darts == any_dm.num_darts


def test_double_dual_is_isomorphic_via_alpha(any_dm) -> None:
    twice = dual(dual(any_dm).dm).dm
    assert is_isomorphism(any_dm, twice, any_dm.alpha)


def test_dual_is_achiral_it_commutes_with_mirror(any_dm) -> None:
    # dual(mirror(M)) = mirror(dual(M)) up to the relabelling alpha.
    left = dual(mirror(any_dm)).dm
    right = mirror(dual(any_dm).dm)
    assert is_isomorphism(left, right, any_dm.alpha)


def test_dual_is_label_invariant(irregular_dm) -> None:
    # Renaming the input darts by pi must rename the output darts by pi.
    # This catches dependence on dart numbering even when M is symmetric.
    rng = random.Random(0)
    perm = list(range(irregular_dm.num_darts))
    rng.shuffle(perm)

    renamed = relabel(irregular_dm, perm)
    assert is_isomorphism(
        dual(irregular_dm).dm, dual(renamed).dm, perm
    )


# ---------------------------------------------------------------------------
# dual: provenance
# ---------------------------------------------------------------------------


def _dart_sets(dm: DartMap, kind: str, ids: list[int]) -> list[frozenset]:
    """Dart sets of the old objects named by canonical ids."""
    orbits = dm.face_orbits() if kind == "face" else dm.vertex_orbits()
    by_min = {min(o): frozenset(o) for o in orbits}
    return [by_min[i] for i in ids]


def test_dual_provenance_vertices_come_from_faces(any_dm) -> None:
    result = dual(any_dm)
    origins = result.vertex_origin

    assert len(origins) == len(result.dm.vertex_orbits())
    assert all(o.kind == "face" for o in origins)
    # canonical id = smallest dart of the old face; each face used once
    assert sorted(o.old_id for o in origins) == sorted(
        min(orbit) for orbit in any_dm.face_orbits()
    )
    # new vertex i really is the old face it names
    for origin, orbit in zip(
        origins, result.dm.vertex_orbits(), strict=True
    ):
        assert _dart_sets(any_dm, "face", [origin.old_id]) == [
            frozenset(orbit)
        ]


def test_dual_provenance_faces_come_from_vertices(any_dm) -> None:
    result = dual(any_dm)
    origins = result.face_origin

    assert len(origins) == len(result.dm.face_orbits())
    assert all(o.kind == "vertex" for o in origins)
    assert sorted(o.old_id for o in origins) == sorted(
        min(orbit) for orbit in any_dm.vertex_orbits()
    )
    # new face j is the set of darts pointing into the old vertex it names
    for origin, orbit in zip(origins, result.dm.face_orbits(), strict=True):
        (old_vertex,) = _dart_sets(any_dm, "vertex", [origin.old_id])
        assert frozenset(any_dm.alpha[d] for d in old_vertex) == frozenset(
            orbit
        )


def test_dual_provenance_darts_are_the_same_darts(any_dm) -> None:
    result = dual(any_dm)
    assert len(result.dart_origin) == any_dm.num_darts
    for d, origin in enumerate(result.dart_origin):
        assert origin == Origin("dart", d, "dual")


def test_dual_merges_nothing(any_dm) -> None:
    # No two raw keys are identified: each new vertex has exactly one.
    result = dual(any_dm)
    assert len(result.merged_from) == len(result.dm.vertex_orbits())
    assert all(len(raw) == 1 for raw in result.merged_from)


def test_dual_provenance_is_label_invariant(irregular_dm) -> None:
    # The *dart sets* named by the provenance move with the relabelling,
    # even though canonical ids (smallest darts) do not.
    rng = random.Random(1)
    perm = list(range(irregular_dm.num_darts))
    rng.shuffle(perm)
    renamed = relabel(irregular_dm, perm)

    def vertex_faces(dm: DartMap) -> set[frozenset]:
        res = dual(dm)
        return set(
            _dart_sets(dm, "face", [o.old_id for o in res.vertex_origin])
        )

    expected = {
        frozenset(perm[d] for d in s) for s in vertex_faces(irregular_dm)
    }
    assert vertex_faces(renamed) == expected


# ---------------------------------------------------------------------------
# dual: orientation (needs coordinates, pure Python)
# ---------------------------------------------------------------------------

Vec = tuple[float, float, float]


def _sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _newell(points: list[Vec]) -> Vec:
    """Area-weighted normal of a polygon, CCW seen from the tip."""
    nx = ny = nz = 0.0
    for i, p in enumerate(points):
        q = points[(i + 1) % len(points)]
        nx += (p[1] - q[1]) * (p[2] + q[2])
        ny += (p[2] - q[2]) * (p[0] + q[0])
        nz += (p[0] - q[0]) * (p[1] + q[1])
    return (nx, ny, nz)


def _affine(points: list[Vec]) -> list[Vec]:
    """Apply a fixed asymmetric invertible affine map (breaks symmetry)."""
    m = ((1.0, 0.3, 0.1), (0.0, 1.7, 0.4), (0.2, 0.0, 2.3))
    t = (0.5, -0.2, 0.9)
    return [
        tuple(_dot(row, p) + t[i] for i, row in enumerate(m))  # type: ignore
        for p in points
    ]


def _cube() -> tuple[list[Vec], list[list[int]]]:
    pts = [
        (float(x), float(y), float(z))
        for z in (0, 1)
        for y in (0, 1)
        for x in (0, 1)
    ]
    faces = [
        [0, 1, 3, 2],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 3, 7, 6],
        [0, 2, 6, 4],
        [1, 3, 7, 5],
    ]
    return pts, faces


def _prism3() -> tuple[list[Vec], list[list[int]]]:
    pts: list[Vec] = [
        (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1),
    ]  # fmt: skip
    faces = [[0, 1, 2], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]]
    return pts, faces


def _pyramid4() -> tuple[list[Vec], list[list[int]]]:
    pts: list[Vec] = [
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.3, 0.4, 1.7),
    ]  # fmt: skip
    faces = [[0, 1, 2, 3], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
    return pts, faces


GEOMETRIC_SOLIDS = [
    ("sheared cube", _cube),
    ("sheared triangular prism", _prism3),
    ("off-centre square pyramid", _pyramid4),
]


def _polar_dual_orientation_signs(
    points: list[Vec], faces: list[list[int]], *, flip: bool
) -> list[float]:
    """Return sign data for every face of the dual of a convex polyhedron.

    The polyhedron is oriented outward-CCW, turned into a ``DartMap``, and
    dualised with :func:`dual`.  Dual vertices are placed by *polarity* about
    an interior point (exact for convex polyhedra), so each dual face must be
    outward-CCW around the old vertex it belongs to.  Returns, per dual face,
    ``dot(newell_normal, old_vertex - centre)``; all positive means the dual
    kept the outward-CCW convention.  ``flip`` mirrors the dual first, as a
    negative control.
    """
    points = _affine(points)
    centre = tuple(sum(p[i] for p in points) / len(points) for i in range(3))

    oriented: list[list[int]] = []
    for face in faces:
        pts = [points[i] for i in face]
        mid = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
        if _dot(_newell(pts), _sub(mid, centre)) < 0:
            face = face[::-1]
        oriented.append(face)

    dm = DartMap.from_face_lists(oriented, num_vertices=len(points))
    # dart k belongs to face dart_face[k] and runs dart_tail[k] -> dart_head[k]
    dart_face: list[int] = []
    dart_tail: list[int] = []
    dart_head: list[int] = []
    for fi, face in enumerate(oriented):
        for i, u in enumerate(face):
            dart_face.append(fi)
            dart_tail.append(u)
            dart_head.append(face[(i + 1) % len(face)])

    polar: list[Vec] = []
    for face in oriented:
        pts = [points[i] for i in face]
        n = _newell(pts)
        mid = tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))
        height = _dot(n, _sub(mid, centre)) / _dot(n, n)
        polar.append(tuple(c / height for c in n))  # type: ignore[arg-type]
        # polar point of the plane {x : n.(x - centre) = |n|^2 * height}

    result = dual(dm).dm
    if flip:
        result = mirror(result)
    signs: list[float] = []
    for orbit in result.face_orbits():
        # dual dart d leaves old face dart_face[d] (sigma* = phi)
        ring = [polar[dart_face[d]] for d in orbit]
        # The old vertex on the left of dual dart d is head(d) = tail(alpha d);
        # in the mirror image left and right swap, so it is tail(d).
        if flip:
            owners = {dart_tail[d] for d in orbit}
        else:
            owners = {dart_tail[dm.alpha[d]] for d in orbit}
        assert len(owners) == 1
        old_vertex = points[owners.pop()]
        signs.append(_dot(_newell(ring), _sub(old_vertex, centre)))
    return signs


@pytest.mark.parametrize(
    ("name", "build"), GEOMETRIC_SOLIDS, ids=[n for n, _ in GEOMETRIC_SOLIDS]
)
def test_dual_preserves_outward_ccw_orientation(name, build) -> None:
    points, faces = build()
    signs = _polar_dual_orientation_signs(points, faces, flip=False)
    assert signs and all(s > 0 for s in signs), name


@pytest.mark.parametrize(
    ("name", "build"), GEOMETRIC_SOLIDS, ids=[n for n, _ in GEOMETRIC_SOLIDS]
)
def test_orientation_test_detects_a_mirrored_dual(name, build) -> None:
    # Negative control: the same check must fail for the mirror image.
    points, faces = build()
    signs = _polar_dual_orientation_signs(points, faces, flip=True)
    assert signs and all(s < 0 for s in signs), name
