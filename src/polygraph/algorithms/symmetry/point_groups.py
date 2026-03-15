"""Point group definitions for polyhedral symmetry classification.

A point group is a finite group of symmetry operations that fix at least one
point.  For convex polyhedra the complete set of relevant families is:

**Polyhedral (exceptional) groups**

* T (order 12) — chiral tetrahedral rotations
* T_d (order 24) — full tetrahedral (T + σ_d mirrors)
* T_h (order 24) — pyritohedral (T + inversion i; no σ_d mirrors)
* O (order 24) — chiral octahedral rotations
* O_h (order 48) — full octahedral (O + inversion)
* I (order 60) — chiral icosahedral rotations
* I_h (order 120) — full icosahedral (I + inversion)

**Dihedral family**

* D_n (order 2n) — n-fold axis + n perpendicular C_2 axes; no mirrors
* D_nh (order 4n) — D_n + horizontal mirror σ_h
* D_nd (order 4n) — D_n + dihedral mirrors σ_d (between the C_2 axes)

**Cyclic family**

* C_n (order n) — pure n-fold rotation axis
* C_nv (order 2n) — C_n + n vertical mirrors σ_v
* C_nh (order 2n) — C_n + horizontal mirror σ_h

**Improper/special groups**

* S_2n (order 2n) — improper rotation axis; S_2 = C_i (inversion)
* C_s (order 2) — single mirror plane; equivalent to C_1h = C_1v
* C_i (order 2) — inversion centre only; equivalent to S_2

Groups are identified by their Schönflies symbol and group order.  No
geometric information (rotation axes, mirror planes) is stored here; the
classification module derives group membership combinatorially from the
dart-map automorphism group.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PointGroup:
    """Named point group with its Schönflies symbol and order.

    Parameters
    ----------
    name : str
        Schönflies notation, e.g. ``"T_d"``, ``"D_3h"``, ``"C_4v"``.
    order : int
        ``|G|``, the number of symmetry operations in the group.
    """

    name: str
    order: int

    def __repr__(self) -> str:
        """Return unambiguous string representation."""
        return f"PointGroup({self.name!r}, order={self.order})"


# ---------------------------------------------------------------------------
# Polyhedral (exceptional) group factories
# ---------------------------------------------------------------------------


def tetrahedral(full: bool = True) -> PointGroup:
    """Return the tetrahedral point group T_d or T.

    Parameters
    ----------
    full : bool
        If ``True`` (default), return T_d (order 24, includes σ_d mirror
        planes and S_4 axes).  If ``False``, return T (order 12, pure
        rotations only — the chiral subgroup).

    Returns
    -------
    PointGroup
        ``T_d`` (order 24) or ``T`` (order 12).

    Notes
    -----
    T_d is the full symmetry group of the regular tetrahedron; its center is
    trivial (no central inversion).  T is the chiral subgroup of T_d and also
    the orientation-preserving subgroup shared by both T_d and T_h.

    See Also
    --------
    pyritohedral : The third tetrahedral group, T_h (order 24, has inversion).
    """
    if full:
        return PointGroup(name="T_d", order=24)
    return PointGroup(name="T", order=12)


def pyritohedral() -> PointGroup:
    """Return the pyritohedral point group T_h (order 24).

    Returns
    -------
    PointGroup
        ``T_h`` with order 24.

    Notes
    -----
    T_h = T × {E, i} — the tetrahedral rotation group augmented by the
    inversion centre i.  Unlike T_d, T_h has an inversion centre (making i
    a central element) but no S_4 improper rotation axes and no σ_d mirrors.

    The orientation-preserving subgroup of T_h is T (order 12), the same as
    for T_d, but the improper elements differ: T_d contains σ_d and S_4
    (neither central), while T_h contains i and S_6 (i is central).

    T_h is the symmetry group of the pyritohedron (an irregular dodecahedron
    with 12 pentagonal faces of T_h symmetry in 3D).

    See Also
    --------
    tetrahedral : T_d (order 24, no inversion) and T (order 12, chiral).
    """
    return PointGroup(name="T_h", order=24)


def octahedral(full: bool = True) -> PointGroup:
    """Return the octahedral point group O_h or O.

    Parameters
    ----------
    full : bool
        If ``True`` (default), return O_h (order 48, includes inversion and
        mirror planes).  If ``False``, return O (order 24, pure rotations
        only — the chiral subgroup).

    Returns
    -------
    PointGroup
        ``O_h`` (order 48) or ``O`` (order 24).

    Notes
    -----
    O_h is the full symmetry group of the cube and the regular octahedron.
    O is its orientation-preserving chiral subgroup.  O arises for the snub
    cube (geometrically), though that solid is not yet implemented.
    """
    if full:
        return PointGroup(name="O_h", order=48)
    return PointGroup(name="O", order=24)


def icosahedral(full: bool = True) -> PointGroup:
    """Return the icosahedral point group I_h or I.

    Parameters
    ----------
    full : bool
        If ``True`` (default), return I_h (order 120, includes inversion and
        mirror planes).  If ``False``, return I (order 60, pure rotations
        only — the chiral subgroup).

    Returns
    -------
    PointGroup
        ``I_h`` (order 120) or ``I`` (order 60).

    Notes
    -----
    I_h is the full symmetry group of the regular dodecahedron and icosahedron.
    I is its chiral subgroup, arising for the snub dodecahedron (not yet
    implemented).
    """
    if full:
        return PointGroup(name="I_h", order=120)
    return PointGroup(name="I", order=60)


# ---------------------------------------------------------------------------
# Dihedral group factories
# ---------------------------------------------------------------------------


def dihedral(n: int, variant: str = "h") -> PointGroup:
    """Return a full dihedral point group D_nh or D_nd (order 4n).

    Parameters
    ----------
    n : int
        Fold order of the principal axis; n ≥ 2.
    variant : str
        ``"h"`` for D_nh (has a horizontal mirror σ_h — prisms, dipyramids),
        ``"d"`` for D_nd (has dihedral mirrors σ_d but no σ_h — antiprisms).

    Returns
    -------
    PointGroup
        ``D_nh`` or ``D_nd`` with order ``4n``.

    Notes
    -----
    D_nh and D_nd are both the full symmetry groups of the respective
    polyhedra and share the chiral dihedral subgroup D_n.  See
    :func:`dihedral_chiral` for the rotation-only variant.
    """
    if n < 2:
        raise ValueError(
            f"n must be >= 2 for dihedral point groups, got {n}"
        )
    if variant not in ("h", "d"):
        raise ValueError(f"variant must be 'h' or 'd', got {variant!r}")
    return PointGroup(name=f"D_{n}{variant}", order=4 * n)


def dihedral_chiral(n: int) -> PointGroup:
    """Return the chiral dihedral group D_n (order 2n).

    Parameters
    ----------
    n : int
        Fold order of the principal axis; n ≥ 2.

    Returns
    -------
    PointGroup
        ``D_n`` with order ``2n``.

    Notes
    -----
    D_n contains one n-fold rotation axis and n perpendicular C_2 axes, but
    no mirror planes and no improper rotations.  It is the orientation-
    preserving (chiral) subgroup shared by both D_nh and D_nd.
    """
    if n < 2:
        raise ValueError(
            f"n must be >= 2 for dihedral point groups, got {n}"
        )
    return PointGroup(name=f"D_{n}", order=2 * n)


# ---------------------------------------------------------------------------
# Cyclic group factories
# ---------------------------------------------------------------------------


def cyclic(n: int) -> PointGroup:
    """Return the C_nv point group (n-fold axis + n vertical mirrors).

    Parameters
    ----------
    n : int
    Fold order of the principal axis; n ≥ 1.

    Returns
    -------
    PointGroup
    ``C_nv`` with order ``2n``.

    Notes
    -----
    C_nv is the symmetry group of a right pyramid with a regular n-gon base:
    one n-fold rotation axis and n vertical mirror planes, giving order 2n.
    """
    if n < 1:
    raise ValueError(f"n must be >= 1 for cyclic groups, got {n}")
    return PointGroup(name=f"C_{n}v", order=2 * n)


def cyclic_pure(n: int) -> PointGroup:
    """Return the cyclic rotation group C_n (order n).

    Parameters
    ----------
    n : int
        Fold order of the rotation axis; n ≥ 1.

    Returns
    -------
    PointGroup
        ``C_n`` with order ``n``.

    Notes
    -----
    C_n contains only n-fold rotations about a single axis; no mirrors, no
    improper rotations.  It is the chiral subgroup of both C_nv and C_nh.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1 for cyclic groups, got {n}")
    return PointGroup(name=f"C_{n}", order=n)


def cyclic_horizontal(n: int) -> PointGroup:
    """Return C_nh (C_n augmented with a horizontal mirror σ_h), order 2n.

    Parameters
    ----------
    n : int
        Fold order of the rotation axis; n ≥ 1.

    Returns
    -------
    PointGroup
        ``C_nh`` with order ``2n``.

    Notes
    -----
    C_nh = C_n × {E, σ_h}.  It contains n rotations and n improper
    rotations (S_n elements).  C_1h = C_s (mirror plane only).
    """
    if n < 1:
        raise ValueError(f"n must be >= 1 for cyclic groups, got {n}")
    return PointGroup(name=f"C_{n}h", order=2 * n)


# ---------------------------------------------------------------------------
# Improper rotation and special groups
# ---------------------------------------------------------------------------


def improper(n: int) -> PointGroup:
    """Return S_2n, an improper rotation group of order 2n.

    Parameters
    ----------
    n : int
        Half the periodicity of the improper axis; n ≥ 1.
        The resulting group is ``S_{2n}`` with order ``2n``.

    Returns
    -------
    PointGroup
        ``S_{2n}`` with order ``2n``.

    Notes
    -----
    S_2n is generated by a single 2n-fold improper rotation (a rotation by
    π/n followed by reflection in a perpendicular plane).  Special cases:
    S_2 = C_i (inversion only), S_4, S_6 (= C_3i), etc.

    Note that ``improper(1)`` returns S_2 = C_i, not S_1 = C_s.
    Use the module-level constants :data:`C_S` and :data:`C_I` for the
    order-2 mirror and inversion groups, respectively.
    """
    if n < 1:
        raise ValueError(
            f"n must be >= 1 for improper rotation groups, got {n}"
        )
    return PointGroup(name=f"S_{2 * n}", order=2 * n)


#: C_s — the mirror-plane group (single σ reflection), order 2.
#: Equivalent to C_1h and C_1v.
C_S = PointGroup(name="C_s", order=2)

#: C_i — the inversion group (single inversion centre i), order 2.
#: Equivalent to S_2.
C_I = PointGroup(name="C_i", order=2)
