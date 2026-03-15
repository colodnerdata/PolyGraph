"""Point group classification for combinatorial dart maps.

Given a dart map and the generators of its automorphism group (from
:mod:`polygraph.algorithms.symmetry.automorphisms`), :func:`classify_symmetry`
determines which named point group the automorphism group belongs to.

Classification uses three combinatorial invariants — no geometry required:

* **Group order** ``|Aut(dm)|`` — narrows candidates to a small family.
* **Vertex orbit count** ``nv`` — prisms have nv=1 (the horizontal mirror
  maps the top vertex ring onto the bottom); dipyramids have nv=2 (apex
  vertices form one orbit, equatorial vertices another).
* **Face orbit count** ``nf`` and **lateral face size** — distinguishes D_nh
  (prism: quad laterals) from D_nd (antiprism: triangular laterals).

See :ref:`classification-table` in the Notes section of
:func:`classify_symmetry` for the full decision table.
"""

from __future__ import annotations

from polygraph.algorithms.symmetry.automorphisms import (
    automorphism_group_order,
    is_orientation_preserving,
)
from polygraph.algorithms.symmetry.orbits import face_orbits, vertex_orbits
from polygraph.algorithms.symmetry.point_groups import (
    PointGroup,
    cyclic,
    dihedral,
    icosahedral,
    octahedral,
    pyritohedral,
    tetrahedral,
)
from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation
from polygraph.structures.traversal import face_darts


class UnknownSymmetryError(ValueError):
    """Raised when the automorphism group cannot be matched to a known group.

    This may occur for polyhedra with symmetry groups not yet covered by
    the classifier (e.g. chiral solids or very unusual topologies).
    """


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _face_sizes_per_orbit(
    generators: list[Permutation], dm: DartMap
) -> list[int]:
    """Return the face size for one representative face per orbit.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for Aut(dm).
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[int]
        One integer per face orbit — the number of boundary darts (= edges =
        vertices) of any face in that orbit.
    """
    return [
        sum(1 for _ in face_darts(dm, orb[0]))
        for orb in face_orbits(generators, dm)
    ]


def _has_reflection(
    generators: list[Permutation], dm: DartMap
) -> bool:
    """Return True if any generator is orientation-reversing.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for Aut(dm).
    dm : DartMap
        Input dart map.

    Returns
    -------
    bool
        ``True`` if at least one generator reverses face orientation.
    """
    return any(not is_orientation_preserving(g, dm) for g in generators)


def _center_has_reflection(
    generators: list[Permutation], dm: DartMap
) -> bool:
    """Return True if the centre contains an orientation-reversing element.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for Aut(dm).
    dm : DartMap
        Input dart map.

    Returns
    -------
    bool
        ``True`` if some orientation-reversing group element commutes with
        every generator (i.e. it lies in the centre of the group).

    Notes
    -----
    This is the combinatorial discriminant between T_d and T_h.  Both groups
    have order 24 and the same orbit structure (flag-transitive), but:

    * T_h = T × {E, i}: the inversion ``i`` is a central element of order 2
      that reverses orientation → ``_center_has_reflection`` returns True.
    * T_d: all orientation-reversing elements (σ_d, S_4) are non-central
      → ``_center_has_reflection`` returns False.

    The check enumerates the full group via BFS (at most 24 elements for the
    cases where this function is called), so performance is not a concern.
    """
    if not generators:
        return False

    n = dm.num_darts
    identity = Permutation.identity(n)
    seen: set[tuple[int, ...]] = {identity.mapping}
    queue: list[Permutation] = [identity]
    all_gens = generators + [g.inverse() for g in generators]

    while queue:
        current = queue.pop()
        for g in all_gens:
            nxt = current.compose(g)
            if nxt.mapping not in seen:
                seen.add(nxt.mapping)
                queue.append(nxt)

    for raw in seen:
        elem = Permutation.from_sequence(list(raw))
        if is_orientation_preserving(elem, dm):
            continue
        # Central iff it commutes with every generator
        if all(
            elem.compose(g).mapping == g.compose(elem).mapping
            for g in generators
        ):
            return True
    return False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def classify_symmetry(
    generators: list[Permutation],
    dm: DartMap,
    order: int | None = None,
) -> PointGroup:
    """Classify the automorphism group of *dm* as a named point group.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for ``Aut(dm)``, as returned by
        :func:`~polygraph.algorithms.symmetry.automorphisms.compute_automorphism_generators`.
    dm : DartMap
        Input dart map.
    order : int, optional
        Pre-computed ``|Aut(dm)|``.  If supplied, the BFS closure inside
        :func:`~polygraph.algorithms.symmetry.automorphisms.automorphism_group_order`
        is skipped (useful when the order is already available from pynauty).

    Returns
    -------
    PointGroup
        Named point group with Schönflies symbol and group order.

    Raises
    ------
    :class:`UnknownSymmetryError`
        If the automorphism group cannot be matched to a supported point group.

    Notes
    -----
    .. _classification-table:

    The classifier uses three combinatorial invariants:

    ======  ===  ===  =================================  ========
    Order   nv   nf   Additional condition               Result
    ======  ===  ===  =================================  ========
    120     1    1    —                                  I_h
    60      1    1    —                                  I
    48      1    1    —                                  O_h
    24      1    1    no reflection                      O
    24      1    1    central reflection present         T_h
    24      1    1    reflection, no central reflection  T_d
    4n      2    1    —                                  D_nh  (dipyramid)
    4n      1    2    lateral face size ≠ 3              D_nh  (prism)
    4n      1    2    lateral face size = 3              D_nd  (antiprism)
    2n      2    2    —                                  C_nv  (pyramid)
    ======  ===  ===  =================================  ========

    *nv* = number of vertex orbits, *nf* = number of face orbits.

    Key orbit-count facts that follow from the symmetry operations:

    * Prisms (D_nh): the horizontal mirror σ_h maps the top ring onto the
      bottom ring, so all 2n vertices are equivalent → nv = 1.
    * Dipyramids (D_nh): no such mapping exists for the apex vertices —
      the two apices form one orbit and the n equatorial vertices another
      → nv = 2.
    * The ``order`` parameter is separated from the Platonic solid tests by
      the (nv, nf) conditions so that, e.g., prism(6) (order = 24, D_6h) and
      the tetrahedron (order = 24, T_d) are not confused.
    """
    if order is None:
        order = automorphism_group_order(generators, dm.num_darts)

    nv = len(vertex_orbits(generators, dm))
    nf = len(face_orbits(generators, dm))

    # --- Icosahedral (dodecahedron, icosahedron) ---
    if order == 120 and nv == 1 and nf == 1:
        return icosahedral(full=True)
    if order == 60 and nv == 1 and nf == 1:
        return icosahedral(full=False)

    # --- Octahedral (cube, octahedron) ---
    if order == 48 and nv == 1 and nf == 1:
        return octahedral(full=True)

    # --- Order-24 flag-transitive: T_d vs T_h vs O ---
    # prism(6) and antiprism(6) also have order=24, but their nf=2, so the
    # nf==1 guard here is necessary and sufficient for this three-way split.
    if order == 24 and nv == 1 and nf == 1:
        if not _has_reflection(generators, dm):
            return octahedral(full=False)           # chiral octahedral → O
        if _center_has_reflection(generators, dm):
            return pyritohedral()                   # pyritohedral → T_h
        return tetrahedral(full=True)               # tetrahedral → T_d

    # --- Dihedral family: order = 4n ---
    if order % 4 == 0:
        n = order // 4
        if n >= 2:
            # Dipyramid: nv=2 (apices ≠ equatorial), nf=1 (all triangles)
            if nv == 2 and nf == 1:
                return dihedral(n, "h")

            # Prism / antiprism branch: nv=1, nf=2 (cap orbit + lateral orbit)
            if nv == 1 and nf == 2:
                face_szs = _face_sizes_per_orbit(generators, dm)
                # n is always the cap face size; filter it out to get lateral.
                # Works for prism(3): n=3, sizes=[3,4] → lateral=4 → D_3h.
                # Works for antiprism(4): n=4, sizes=[4,3] → lateral=3 → D_4d.
                laterals = [s for s in face_szs if s != n]
                if laterals:
                    lateral_size = laterals[0]
                    if lateral_size == 3:
                        return dihedral(n, "d")   # antiprism → D_nd
                    return dihedral(n, "h")       # prism → D_nh

    # --- Cyclic-pyramidal: order = 2n ---
    if order % 2 == 0 and nv == 2 and nf == 2:
        n = order // 2
        return cyclic(n)   # pyramid → C_nv

    raise UnknownSymmetryError(
        f"Cannot classify automorphism group: "
        f"|Aut|={order}, vertex_orbits={nv}, face_orbits={nf}"
    )
