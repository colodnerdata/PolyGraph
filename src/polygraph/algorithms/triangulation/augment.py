"""Barycentric subdivision of a combinatorial dart map.

The barycentric subdivision (order complex of the face poset) produces a
canonical, coordinate-free triangulation that automatically respects every
automorphism of the original map.  Each original dart ``d`` generates two
*flags* — (vertex, edge, face) triples — yielding one triangle per flag.
The construction is uniform: every face, including triangles, is subdivided.

New dart indexing
-----------------
For each original dart ``d`` in ``0..n-1``, six new darts ``6d..6d+5`` are
created:

* ``6d+0, 6d+1, 6d+2`` — boundary of the *start-flag* triangle
  (vertex → edge-midpoint → face-center → vertex)
* ``6d+3, 6d+4, 6d+5`` — boundary of the *end-flag* triangle
  (far-vertex → face-center → edge-midpoint → far-vertex)

Counts (original map has V vertices, E edges, F faces, n = 2E darts):

* Vertices:  V + E + F
* Edges:     6E
* Faces:     4E  (all triangles)
* Darts:     12E (= 6n)
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto

from polygraph.structures.dart_map import DartMap, _inverse_permutation


class CellType(Enum):
    """Type of original cell a barycentric vertex came from."""

    VERTEX = auto()
    EDGE = auto()
    FACE = auto()


@dataclass(frozen=True)
class CellOrigin:
    """Maps a new vertex back to the original cell it represents.

    Attributes
    ----------
    cell_type : CellType
        Whether this vertex represents an original vertex, edge-midpoint,
        or face-center.
    representative : int
        A representative dart of the original cell (vertex-rep, edge-rep,
        or face-rep following the dart-map conventions).
    """

    cell_type: CellType
    representative: int


@dataclass(frozen=True)
class BarycentricResult:
    """Result of a barycentric subdivision.

    Attributes
    ----------
    dart_map : DartMap
        The subdivided, all-triangular dart map.
    cell_map : dict[int, CellOrigin]
        Maps each new vertex representative dart to the original cell it
        came from.
    """

    dart_map: DartMap
    cell_map: dict[int, CellOrigin]


def barycentric_subdivision(dm: DartMap) -> BarycentricResult:
    """Compute the barycentric subdivision of a dart map.

    Parameters
    ----------
    dm : DartMap
        Input combinatorial map (any genus, any face sizes).

    Returns
    -------
    BarycentricResult
        The subdivided dart map (all faces are triangles) together with a
        cell map recording the origin of each new vertex.

    Notes
    -----
    The subdivision is canonical — it uses only the combinatorial
    permutations sigma, alpha, phi — so every automorphism of ``dm``
    extends to an automorphism of the result.
    """
    n = dm.num_darts
    new_n = 6 * n

    # Pre-compute phi and phi_inv as arrays for speed.
    alpha = dm.alpha
    sigma_inv = _inverse_permutation(sigma)
    phi = [sigma_inv[alpha[d]] for d in range(n)]
    phi_inv = _inverse_permutation(phi)

    # ------------------------------------------------------------------
    # Build new_alpha  (edge-pairing involution, size 6n)
    # ------------------------------------------------------------------
    new_alpha = [0] * new_n
    for d in range(n):
        base = 6 * d

        # Type 1: e–f edge (between start & end flag of same dart d)
        new_alpha[base + 1] = base + 4
        new_alpha[base + 4] = base + 1

        # Type 2: v–e edge (change face, cross original edge)
        a_d = alpha[d]
        new_alpha[base + 0] = 6 * a_d + 5
        new_alpha[base + 5] = 6 * a_d + 0

        # Type 3: v–f edge (change edge, stay in same face)
        new_alpha[base + 2] = 6 * phi_inv[d] + 3
        new_alpha[base + 3] = 6 * phi[d] + 2

    # ------------------------------------------------------------------
    # Build new_phi  (face boundary, size 6n — each triangle is a 3-cycle)
    # ------------------------------------------------------------------
    new_phi = [0] * new_n
    for d in range(n):
        base = 6 * d
        # Start-flag triangle: 0 → 1 → 2 → 0
        new_phi[base + 0] = base + 1
        new_phi[base + 1] = base + 2
        new_phi[base + 2] = base + 0
        # End-flag triangle: 3 → 4 → 5 → 3
        new_phi[base + 3] = base + 4
        new_phi[base + 4] = base + 5
        new_phi[base + 5] = base + 3

    # ------------------------------------------------------------------
    # Derive new_sigma = new_alpha ∘ new_phi⁻¹
    # ------------------------------------------------------------------
    new_phi_inv = _inverse_permutation(new_phi)
    new_sigma = [new_alpha[new_phi_inv[i]] for i in range(new_n)]

    subdivided = DartMap(sigma=new_sigma, alpha=new_alpha)

    # ------------------------------------------------------------------
    # Build cell_map: new vertex representative → original cell
    # ------------------------------------------------------------------
    cell_map = _build_cell_map(dm, subdivided)

    return BarycentricResult(dart_map=subdivided, cell_map=cell_map)


def _build_cell_map(
    original: DartMap, subdivided: DartMap
) -> dict[int, CellOrigin]:
    """Identify which original cell each new vertex corresponds to.

    New vertex types by new-dart offset within the 6-dart block of
    original dart ``d``:

    * ``6d+0`` starts at the original vertex of ``d``  →  VERTEX
    * ``6d+1`` starts at the edge-midpoint of ``d``    →  EDGE
    * ``6d+2`` starts at the face-center of ``d``      →  FACE
    * ``6d+3`` starts at the far vertex of ``d``        →  VERTEX
    * ``6d+4`` starts at the face-center of ``d``      →  FACE
    * ``6d+5`` starts at the edge-midpoint of ``d``    →  EDGE
    """

    alpha = original.alpha

    # Precompute canonical representatives for original vertices and faces:
    # - Vertex rep: min dart in the original sigma orbit.
    # - Face rep:   min dart in the original phi orbit.
    vertex_rep: list[int] = list(range(n_orig))
    for orbit in _cycle_orbits(original.sigma):
        rep = min(orbit)
        for d in orbit:
            vertex_rep[d] = rep

    face_rep: list[int] = list(range(n_orig))
    for orbit in _cycle_orbits(original.phi):
        rep = min(orbit)
        for d in orbit:
            face_rep[d] = rep

    # For each new dart, determine (cell_type, original_representative).
    def _cell_of_new_dart(new_d: int) -> CellOrigin:
        d = new_d // 6
        offset = new_d % 6
        if offset == 0:
            # original vertex of d — use canonical vertex rep
            return CellOrigin(CellType.VERTEX, vertex_rep[d])
        elif offset == 1:
            # edge-midpoint — use canonical edge rep min(d, alpha[d])
            return CellOrigin(CellType.EDGE, min(d, alpha[d]))
        elif offset == 2:
            # face-center — use canonical face rep
            return CellOrigin(CellType.FACE, face_rep[d])
        elif offset == 3:
            # far vertex = vertex of alpha(d) — use its canonical vertex rep
            return CellOrigin(CellType.VERTEX, vertex_rep[alpha[d]])
        elif offset == 4:
            # face-center (same face as d) — use canonical face rep
            return CellOrigin(CellType.FACE, face_rep[d])
        else:  # offset == 5
            # edge-midpoint (same edge as d)
            return CellOrigin(CellType.EDGE, min(d, alpha[d]))

    # Walk sigma orbits of the *subdivided* map to find vertex reps,
    # then annotate each rep with its original cell.
    new_sigma = subdivided.sigma
    seen = [False] * subdivided.num_darts
    cell_map: dict[int, CellOrigin] = {}
    for start in range(subdivided.num_darts):
        if seen[start]:
            continue
        # start is the representative of this new vertex
        cell_map[start] = _cell_of_new_dart(start)
        cur = start
        while True:
            seen[cur] = True
            cur = new_sigma[cur]
            if cur == start:
                break

    return cell_map
