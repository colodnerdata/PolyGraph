"""Orbit computations for permutation groups acting on DartMap entities."""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation


def compute_orbits(generators: list[Permutation], n: int) -> list[list[int]]:
    """Compute element orbits under a generated permutation group.

    Parameters
    ----------
    generators : list[Permutation]
        Group generators acting on ``0..n-1``.
    n : int
        Domain size.

    Returns
    -------
    list[list[int]]
        Partition of ``0..n-1`` into sorted orbits.

    Raises
    ------
    ValueError
        If any generator does not act on ``n`` elements.
    """
    for generator in generators:
        if len(generator) != n:
            raise ValueError("all generators must have size n")

    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        root_a = find(a)
        root_b = find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    for generator in generators:
        for x in range(n):
            union(x, generator[x])

    parts: dict[int, list[int]] = {}
    for x in range(n):
        parts.setdefault(find(x), []).append(x)

    return [sorted(orbit) for _, orbit in sorted(parts.items())]


def dart_orbits(generators: list[Permutation], dm: DartMap) -> list[list[int]]:
    """Compute dart orbits under ``generators``."""
    return compute_orbits(generators, dm.num_darts)


def vertex_orbits(generators: list[Permutation], dm: DartMap) -> list[list[int]]:
    """Compute vertex-orbit partition induced by dart automorphisms."""
    vertex_reps = [min(orbit) for orbit in dm.vertex_orbits()]
    rep_of_dart: dict[int, int] = {}
    for orbit in dm.vertex_orbits():
        rep = min(orbit)
        for dart in orbit:
            rep_of_dart[dart] = rep

    index_of_rep = {rep: i for i, rep in enumerate(vertex_reps)}
    action: list[Permutation] = []
    for generator in generators:
        action.append(
            Permutation.from_sequence(
                [index_of_rep[rep_of_dart[generator[rep]]] for rep in vertex_reps]
            )
        )
    rep_orbits = compute_orbits(action, len(vertex_reps))
    return [[vertex_reps[i] for i in orbit] for orbit in rep_orbits]


def edge_orbits(generators: list[Permutation], dm: DartMap) -> list[list[int]]:
    """Compute edge-orbit partition induced by dart automorphisms."""
    edge_reps = sorted({min(d, dm.alpha[d]) for d in range(dm.num_darts)})
    index_of_rep = {rep: i for i, rep in enumerate(edge_reps)}

    action: list[Permutation] = []
    for generator in generators:
        action.append(
            Permutation.from_sequence(
                [
                    index_of_rep[min(generator[rep], dm.alpha[generator[rep]])]
                    for rep in edge_reps
                ]
            )
        )
    rep_orbits = compute_orbits(action, len(edge_reps))
    return [[edge_reps[i] for i in orbit] for orbit in rep_orbits]


def face_orbits(generators: list[Permutation], dm: DartMap) -> list[list[int]]:
    """Compute face-orbit partition induced by dart automorphisms."""
    face_reps = [min(orbit) for orbit in dm.face_orbits()]
    rep_of_dart: dict[int, int] = {}
    for orbit in dm.face_orbits():
        rep = min(orbit)
        for dart in orbit:
            rep_of_dart[dart] = rep

    index_of_rep = {rep: i for i, rep in enumerate(face_reps)}
    action: list[Permutation] = []
    for generator in generators:
        action.append(
            Permutation.from_sequence(
                [index_of_rep[rep_of_dart[generator[rep]]] for rep in face_reps]
            )
        )
    rep_orbits = compute_orbits(action, len(face_reps))
    return [[face_reps[i] for i in orbit] for orbit in rep_orbits]


__all__ = [
    "compute_orbits",
    "dart_orbits",
    "edge_orbits",
    "face_orbits",
    "vertex_orbits",
]
