"""Core combinatorial map structure based on darts and permutations.

A DartMap encodes the topology of a graph embedding through three fundamental
permutations on the set of darts:
sigma — cyclic order of darts around each vertex
alpha — pairing of opposite darts along each edge (an involution)
phi(d) = sigma^{-1}(alpha(d)) — traversal of darts around faces
This representation cleanly separates combinatorial topology from geometry.
The module provides utilities to construct maps, iterate vertex/edge/face
orbits, and compute topological invariants such as the Euler characteristic
and genus.
The structure is designed to support planar graph algorithms, symmetry
analysis, and polyhedral realizations while remaining independent of any
coordinate system.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field


def _inverse_permutation(perm: Sequence[int]) -> list[int]:
    """Compute the inverse of a permutation.

    Parameters
    ----------
    perm : Sequence[int]
        Permutation in one-line notation, where ``perm[i]`` is the image of
        index ``i``.

    Returns
    -------
    list[int]
        Inverse permutation ``inv`` satisfying ``inv[perm[i]] == i``.

    Raises
    ------
    ValueError
        If ``perm`` is not a valid permutation on ``range(len(perm))``.
    """
    n = len(perm)
    inv = [-1] * n
    for i, j in enumerate(perm):
        if not 0 <= j < n:
            raise ValueError(
                f"Permutation target out of range at index {i}: {j}."
            )
        if inv[j] != -1:
            raise ValueError(
                f"Permutation is not bijective; duplicate target {j}."
            )
        inv[j] = i
    return inv


def _cycle_orbits(perm: Sequence[int]) -> list[list[int]]:
    """Decompose a permutation into disjoint cycles.

    Parameters
    ----------
    perm : Sequence[int]
        Permutation in one-line notation.

    Returns
    -------
    list[list[int]]
        Disjoint cycles of ``perm`` as lists of dart indices.
    """
    n = len(perm)
    seen = [False] * n
    cycles: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        cycle: list[int] = []
        d = start
        while not seen[d]:
            seen[d] = True
            cycle.append(d)
            d = perm[d]
        cycles.append(cycle)
    return cycles


@dataclass(slots=True)
class DartMap:
    """Combinatorial map represented by vertex and edge permutations.

    Parameters
    ----------
    sigma : list[int]
        Vertex rotation permutation where ``sigma[d]`` is the next dart around
        the same vertex.
    alpha : list[int]
        Edge involution where ``alpha[d]`` is the opposite dart along the same
        edge.

    Notes
    -----
    A dart is an oriented half-edge. The face permutation is derived as
    ``phi = sigma^{-1} o alpha``.
    """

    sigma: list[int]
    alpha: list[int]
    _sigma_inv: list[int] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Validate and normalize permutation data.

        Runs after dataclass initialization.

        Raises
        ------
        ValueError
            If permutation sizes mismatch, permutations are invalid, the number
            of darts is odd, or ``alpha`` is not a fixed-point-free involution.
        """
        self.sigma = list(self.sigma)
        self.alpha = list(self.alpha)

        n = len(self.sigma)
        if n == 0:
            raise ValueError("Dart map must contain at least one dart.")
        if len(self.alpha) != n:
            raise ValueError(
                "sigma and alpha must have the same size "
                f"(got {n} and {len(self.alpha)})."
            )

        self._sigma_inv = _inverse_permutation(self.sigma)
        _inverse_permutation(self.alpha)

        if n % 2 != 0:
            raise ValueError(
                "Number of darts must be even for a fixed-point-free alpha."
            )

        for d, twin in enumerate(self.alpha):
            if twin == d:
                raise ValueError(
                    f"alpha must be fixed-point-free; alpha({d}) == {d}."
                )
            if self.alpha[twin] != d:
                raise ValueError(
                    f"alpha must be an involution; alpha(alpha({d})) = "
                    f"{self.alpha[twin]}."
                )

    @property
    def num_darts(self) -> int:
        """Return the number of darts in the map.

        Returns
        -------
        int
            Total number of darts.
        """
        return len(self.sigma)

    @property
    def num_edges(self) -> int:
        """Return the number of undirected edges.

        Returns
        -------
        int
            Number of edges, computed as ``num_darts // 2``.
        """
        return self.num_darts // 2

    def validate_dart(self, d: int) -> None:
        """Validate that a dart index is in range ``[0, num_darts)``.

        Parameters
        ----------
        d : int
            Dart index to validate.

        Raises
        ------
        IndexError
            If ``d`` is not in ``[0, num_darts)``.
        """
        if not 0 <= d < self.num_darts:
            raise IndexError(
                f"Dart index {d} out of range [0, {self.num_darts})."
            )

    def _check_dart(self, d: int) -> None:
        """Validate that a dart index is in range.

        Parameters
        ----------
        d : int
            Dart index to validate.

        Raises
        ------
        IndexError
            If ``d`` is not in ``[0, num_darts)``.
        """
        self.validate_dart(d)

    def phi(self, d: int) -> int:
        """Return the next dart around the incident face.

        Parameters
        ----------
        d : int
            Input dart index.

        Returns
        -------
        int
            Dart index ``phi(d)``.

        Raises
        ------
        IndexError
            If ``d`` is out of range.
        """
        self._check_dart(d)
        return self._sigma_inv[self.alpha[d]]

    def edge_of_dart(self, d: int) -> int:
        """Return a canonical edge identifier for a dart.

        Parameters
        ----------
        d : int
            Input dart index.

        Returns
        -------
        int
            Canonical identifier of the undirected edge containing ``d``,
            represented by ``min(d, alpha[d])``.

        Raises
        ------
        IndexError
            If ``d`` is out of range.
        """
        self._check_dart(d)
        return min(d, self.alpha[d])

    def vertex_orbits(self) -> list[list[int]]:
        """Return vertex cycles induced by ``sigma``.

        Returns
        -------
        list[list[int]]
            Darts grouped into vertex orbits.
        """
        return _cycle_orbits(self.sigma)

    def face_orbits(self) -> list[list[int]]:
        """Return face cycles induced by ``phi``.

        Returns
        -------
        list[list[int]]
            Darts grouped into face orbits.
        """
        phi_perm = [self.phi(d) for d in range(self.num_darts)]
        return _cycle_orbits(phi_perm)

    def connected_components(self) -> list[list[int]]:
        """Return connected components in the dart adjacency graph.

        Returns
        -------
        list[list[int]]
            Components of darts under adjacency generated by ``sigma``,
            ``sigma^{-1}``, and ``alpha``.
        """
        seen = [False] * self.num_darts
        components: list[list[int]] = []
        for start in range(self.num_darts):
            if seen[start]:
                continue
            stack = [start]
            component: list[int] = []
            seen[start] = True
            while stack:
                d = stack.pop()
                component.append(d)
                for nxt in (self.sigma[d], self._sigma_inv[d], self.alpha[d]):
                    if not seen[nxt]:
                        seen[nxt] = True
                        stack.append(nxt)
            components.append(component)
        return components

    def euler_characteristic(self) -> int:
        """Compute the Euler characteristic.

        Returns
        -------
        int
            Euler characteristic ``V - E + F`` for the map.
        """
        vertices = len(self.vertex_orbits())
        edges = self.num_edges
        faces = len(self.face_orbits())
        return vertices - edges + faces

    def genus(self) -> int:
        """Compute orientable genus for a closed surface map.

        Returns
        -------
        int
            Orientable genus ``g`` satisfying ``chi = 2c - 2g``, where ``c`` is
            the number of connected components.

        Raises
        ------
        ValueError
            If the computed genus is non-integral.
        """
        chi = self.euler_characteristic()
        components = len(self.connected_components())
        numerator = 2 * components - chi
        if numerator % 2 != 0:
            raise ValueError(
                "Genus is not integral; map may not represent an orientable "
                "closed surface."
            )
        return numerator // 2

    @classmethod
    def from_face_lists(
        cls, faces: Sequence[Sequence[int]], num_vertices: int
    ) -> DartMap:
        """Construct a map from oriented face vertex cycles.

        Parameters
        ----------
        faces : Sequence[Sequence[int]]
            Faces as cyclic sequences of vertex indices.
        num_vertices : int
            Total number of vertices available for indexing.

        Returns
        -------
        DartMap
            Constructed combinatorial map.
        """
        return dart_map_from_face_lists(faces=faces, num_vertices=num_vertices)


def dart_map_from_face_lists(
    faces: Sequence[Sequence[int]], num_vertices: int
) -> DartMap:
    """Construct a `DartMap` from oriented face boundary vertex lists.

    Parameters
    ----------
    faces : Sequence[Sequence[int]]
        Face boundaries as cyclic vertex index lists.
    num_vertices : int
        Number of vertices in the input indexing scheme.

    Returns
    -------
    DartMap
        Map whose darts correspond to directed face boundary edges.

    Raises
    ------
    ValueError
        If faces are empty, contain invalid vertices or degenerate edges, if
        directed edges repeat, or if opposite directions are missing.

    Notes
    -----
    ``faces`` is expected to describe a closed 2-manifold mesh:
    each undirected edge must appear exactly twice, once in each direction.
    """
    if num_vertices <= 0:
        raise ValueError(f"num_vertices must be positive, got {num_vertices}.")
    if not faces:
        raise ValueError("At least one face is required.")

    directed_to_dart: dict[tuple[int, int], int] = {}
    face_dart_cycles: list[list[int]] = []

    for face_index, face in enumerate(faces):
        face_vertices = list(face)
        if len(face_vertices) < 3:
            raise ValueError(
                f"Face {face_index} has {len(face_vertices)} "
                "vertices; expected >= 3."
            )

        for vertex in face_vertices:
            if not 0 <= vertex < num_vertices:
                raise ValueError(
                    f"Face {face_index} references out-of-range "
                    f"vertex {vertex} (num_vertices={num_vertices})."
                )

        cycle: list[int] = []
        face_len = len(face_vertices)
        for i, u in enumerate(face_vertices):
            v = face_vertices[(i + 1) % face_len]
            if u == v:
                raise ValueError(
                    f"Face {face_index} contains a zero-length "
                    f"edge ({u}, {v})."
                )
            key = (u, v)
            if key in directed_to_dart:
                raise ValueError(
                    f"Directed edge {key} appears more than once; "
                    "faces must define a valid embedding."
                )
            dart = len(directed_to_dart)
            directed_to_dart[key] = dart
            cycle.append(dart)

        face_dart_cycles.append(cycle)

    num_darts = len(directed_to_dart)
    phi = [-1] * num_darts
    for cycle in face_dart_cycles:
        cycle_len = len(cycle)
        for i, d in enumerate(cycle):
            phi[d] = cycle[(i + 1) % cycle_len]

    alpha = [-1] * num_darts
    for (u, v), d in directed_to_dart.items():
        twin = directed_to_dart.get((v, u))
        if twin is None:
            raise ValueError(
                f"Missing opposite directed edge for ({u}, {v}); "
                "mesh is not closed."
            )
        alpha[d] = twin

    phi_inv = _inverse_permutation(phi)
    sigma = [alpha[phi_inv[d]] for d in range(num_darts)]
    return DartMap(sigma=sigma, alpha=alpha)


def edge_of_dart(dm: DartMap, d: int) -> int:
    """Return the canonical edge identifier of a dart.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    d : int
        Input dart index.

    Returns
    -------
    int
        Canonical edge identifier for ``d``.
    """
    return dm.edge_of_dart(d)


def vertex_orbits(dm: DartMap) -> list[list[int]]:
    """Return vertex orbits of a map.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[list[int]]
        Darts grouped by vertex cycles.
    """
    return dm.vertex_orbits()


def face_orbits(dm: DartMap) -> list[list[int]]:
    """Return face orbits of a map.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[list[int]]
        Darts grouped by face cycles.
    """
    return dm.face_orbits()


def euler_characteristic(dm: DartMap) -> int:
    """Return Euler characteristic of a map.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    int
        Euler characteristic ``V - E + F``.
    """
    return dm.euler_characteristic()


def genus(dm: DartMap) -> int:
    """Return orientable genus of a map.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    int
        Orientable genus.
    """
    return dm.genus()
