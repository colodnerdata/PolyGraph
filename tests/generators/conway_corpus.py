"""Irregular convex genus-0 polyhedra and map helpers for operator tests.

Regular solids are a weak test set for combinatorial operators: their large
symmetry groups let a wrong dart pairing still produce a map isomorphic to
the right answer, and every vertex has the same degree and every face the
same size.  This module builds seeded, mostly asymmetric polyhedra using only
the standard library (CI does not install NumPy or SciPy).

By Steinitz's theorem the combinatorial types of convex polyhedra are exactly
the 3-connected planar graphs, so the corpus only has to be 3-connected and
planar.  Construction:

1. Start from the tetrahedron and *stack* a new vertex into a random
   triangular face (planar, 3-connected, degrees become uneven).
2. Apply random *edge flips* to the triangulation (still a simple
   triangulation, so still 3-connected) to destroy stacking structure.
3. Optionally *delete* random edges, merging the two faces they separate,
   whenever the result is still 3-connected.  This yields polyhedra with
   mixed face sizes as well as mixed vertex degrees.

All faces are listed in the same cyclic direction (the orientation convention
of :meth:`DartMap.from_face_lists`), so the result is a valid closed map.
"""

from __future__ import annotations

import random
from collections import deque
from itertools import combinations

from polygraph.structures.dart_map import DartMap

Faces = list[list[int]]

_TETRAHEDRON: Faces = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]]


def relabel(dm: DartMap, perm: list[int]) -> DartMap:
    """Return ``dm`` with dart ``d`` renamed to ``perm[d]``.

    The result ``dm'`` satisfies ``is_isomorphism(dm, dm', perm)``.
    """
    n = dm.num_darts
    sigma = [0] * n
    alpha = [0] * n
    for d in range(n):
        sigma[perm[d]] = perm[dm.sigma[d]]
        alpha[perm[d]] = perm[dm.alpha[d]]
    return DartMap(sigma=sigma, alpha=alpha)


def mirror(dm: DartMap) -> DartMap:
    """Return the mirror image of ``dm`` (``sigma`` becomes its inverse)."""
    inverse = [0] * dm.num_darts
    for d, nxt in enumerate(dm.sigma):
        inverse[nxt] = d
    return DartMap(sigma=inverse, alpha=dm.alpha)


def orientation_preserving_aut_order(dm: DartMap) -> int:
    """Return the number of orientation-preserving automorphisms.

    For a connected map an automorphism is determined by the image of a
    single dart, so it suffices to count the darts ``r`` whose canonical
    breadth-first relabelling from ``r`` equals that from dart ``0``.
    """
    n = dm.num_darts

    def code(root: int) -> tuple[tuple[int, int], ...]:
        label = {root: 0}
        order = [root]
        queue = deque([root])
        while queue:
            d = queue.popleft()
            for nxt in (dm.sigma[d], dm.alpha[d]):
                if nxt not in label:
                    label[nxt] = len(label)
                    order.append(nxt)
                    queue.append(nxt)
        return tuple(
            (label[dm.sigma[d]], label[dm.alpha[d]]) for d in order
        )

    base = code(0)
    return sum(1 for root in range(n) if code(root) == base)


def _directed_edges(faces: Faces) -> dict[tuple[int, int], int]:
    """Map each directed edge ``(u, v)`` to the index of its face."""
    edges: dict[tuple[int, int], int] = {}
    for index, face in enumerate(faces):
        for i, u in enumerate(face):
            edges[(u, face[(i + 1) % len(face)])] = index
    return edges


def _degrees(edges: dict[tuple[int, int], int], n: int) -> list[int]:
    degree = [0] * n
    for u, _ in edges:
        degree[u] += 1
    return degree


def _is_three_connected(faces: Faces, n: int) -> bool:
    """Return whether the graph of ``faces`` is 3-connected (brute force)."""
    adjacency: list[set[int]] = [set() for _ in range(n)]
    for u, v in _directed_edges(faces):
        adjacency[u].add(v)

    def connected_without(removed: tuple[int, ...]) -> bool:
        rest = [v for v in range(n) if v not in removed]
        seen = {rest[0]}
        stack = [rest[0]]
        while stack:
            for w in adjacency[stack.pop()]:
                if w not in removed and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return len(seen) == len(rest)

    if n < 4 or min(len(a) for a in adjacency) < 3:
        return False
    singles = ((v,) for v in range(n))
    return all(
        connected_without(cut)
        for cut in (*singles, *combinations(range(n), 2))
    )


def _stack(faces: Faces, n: int, rng: random.Random) -> int:
    """Insert a vertex into a random triangle; return the new vertex count."""
    index = rng.randrange(len(faces))
    a, b, c = faces[index]
    faces[index] = [a, b, n]
    faces.extend([[b, c, n], [c, a, n]])
    return n + 1


def _try_flip(faces: Faces, n: int, rng: random.Random) -> bool:
    """Flip a random edge of a triangulation if it stays simple."""
    edges = _directed_edges(faces)
    degree = _degrees(edges, n)
    a, b = rng.choice(sorted(edges))
    i1, i2 = edges[(a, b)], edges[(b, a)]
    c = next(v for v in faces[i1] if v not in (a, b))
    d = next(v for v in faces[i2] if v not in (a, b))
    if degree[a] < 4 or degree[b] < 4 or (c, d) in edges or c == d:
        return False
    faces[i1] = [a, d, c]
    faces[i2] = [d, b, c]
    return True


def _try_delete(faces: Faces, n: int, rng: random.Random) -> bool:
    """Delete a random edge, merging its two faces, if still 3-connected."""
    edges = _directed_edges(faces)
    a, b = rng.choice(sorted(edges))
    i1, i2 = edges[(a, b)], edges[(b, a)]

    def rotated(face: list[int], first: int, second: int) -> list[int]:
        for i, u in enumerate(face):
            if u == first and face[(i + 1) % len(face)] == second:
                return face[i:] + face[:i]
        raise AssertionError("edge not on face")

    f1 = rotated(faces[i1], a, b)  # [a, b, p1, ..., pk]
    f2 = rotated(faces[i2], b, a)  # [b, a, q1, ..., qm]
    merged = [a, *f2[2:], b, *f1[2:]]
    candidate = [f for i, f in enumerate(faces) if i not in (i1, i2)]
    candidate.append(merged)
    if not _is_three_connected(candidate, n):
        return False
    faces[:] = candidate
    return True


def irregular_polyhedron(
    seed: int, n_vertices: int = 9, deletions: int = 0
) -> DartMap:
    """Build a seeded irregular convex genus-0 polyhedron.

    Parameters
    ----------
    seed : int
        Seed for the private random generator.
    n_vertices : int
        Number of vertices (at least 4); the edge deletions do not change it.
    deletions : int
        Number of edge-deletion attempts after the triangulation is built;
        each successful one merges two faces.

    Returns
    -------
    DartMap
        A 3-connected planar map, i.e. the combinatorial type of a convex
        polyhedron.
    """
    rng = random.Random(seed)
    faces = [list(face) for face in _TETRAHEDRON]
    n = 4
    while n < n_vertices:
        n = _stack(faces, n, rng)
    for _ in range(8 * n):
        _try_flip(faces, n, rng)
    for _ in range(4 * deletions):
        _try_delete(faces, n, rng)
    return DartMap.from_face_lists(faces, num_vertices=n)


def irregular_corpus() -> list[tuple[str, DartMap]]:
    """Return a deterministic corpus of irregular polyhedra.

    Returns
    -------
    list[tuple[str, DartMap]]
        ``(name, map)`` pairs covering simplicial polyhedra (all faces
        triangles) and polyhedra with mixed face sizes and vertex degrees.
    """
    corpus: list[tuple[str, DartMap]] = []
    for seed in range(14):
        n_vertices = 6 + seed % 7
        deletions = 0 if seed < 5 else 1 + seed % 4
        name = f"irregular(seed={seed}, n={n_vertices}, del={deletions})"
        dm = irregular_polyhedron(seed, n_vertices, deletions)
        corpus.append((name, dm))
    return corpus
