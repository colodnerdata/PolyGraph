"""Point group definitions for polyhedral symmetry classification.

A point group is a finite group of symmetry operations that fix at least one
point.  For convex polyhedra the relevant families are:

* **Polyhedral** — T_d (tetrahedral), O_h (octahedral), I_h (icosahedral)
  and their chiral subgroups T, O, I.
* **Dihedral** — D_nh (prisms, dipyramids) and D_nd (antiprisms).
* **Cyclic-pyramidal** — C_nv (pyramids).

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
# Factory functions
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
    return PointGroup(name=f"C_{n}v", order=2 * n)


def dihedral(n: int, variant: str = "h") -> PointGroup:
    """Return a dihedral point group D_nh or D_nd.

    Parameters
    ----------
    n : int
        Fold order of the principal axis; n ≥ 2.
    variant : str
        ``"h"`` for D_nh (has a horizontal mirror σ_h), ``"d"`` for D_nd
        (has dihedral mirror planes but no σ_h).

    Returns
    -------
    PointGroup
        ``D_nh`` or ``D_nd`` with order ``4n``.

    Notes
    -----
    D_nh is the symmetry group of a right prism or bipyramid with a regular
    n-gon base.  D_nd is the symmetry group of a right antiprism.  Both have
    order 4n, but they differ in whether a horizontal mirror plane exists.
    """
    if variant not in ("h", "d"):
        raise ValueError(f"variant must be 'h' or 'd', got {variant!r}")
    return PointGroup(name=f"D_{n}{variant}", order=4 * n)


def tetrahedral(full: bool = True) -> PointGroup:
    """Return the tetrahedral point group T_d or T.

    Parameters
    ----------
    full : bool
        If ``True`` (default), return T_d (order 24, includes improper
        rotations).  If ``False``, return T (order 12, pure rotations only).

    Returns
    -------
    PointGroup
        ``T_d`` (order 24) or ``T`` (order 12).

    Notes
    -----
    T_d is the full symmetry group of the regular tetrahedron.  T is its
    chiral (orientation-preserving) subgroup of pure rotations.
    """
    if full:
        return PointGroup(name="T_d", order=24)
    return PointGroup(name="T", order=12)


def octahedral(full: bool = True) -> PointGroup:
    """Return the octahedral point group O_h or O.

    Parameters
    ----------
    full : bool
        If ``True`` (default), return O_h (order 48, includes improper
        rotations).  If ``False``, return O (order 24, pure rotations only).

    Returns
    -------
    PointGroup
        ``O_h`` (order 48) or ``O`` (order 24).

    Notes
    -----
    O_h is the full symmetry group of the cube and octahedron.  T is its
    chiral subgroup.
    """
    if full:
        return PointGroup(name="O_h", order=48)
    return PointGroup(name="O", order=24)


def icosahedral(full: bool = True) -> PointGroup:
    """Return the icosahedral point group I_h or I.

    Parameters
    ----------
    full : bool
        If ``True`` (default), return I_h (order 120, includes improper
        rotations).  If ``False``, return I (order 60, pure rotations only).

    Returns
    -------
    PointGroup
        ``I_h`` (order 120) or ``I`` (order 60).

    Notes
    -----
    I_h is the full symmetry group of the dodecahedron and icosahedron.  I is
    its chiral subgroup.
    """
    if full:
        return PointGroup(name="I_h", order=120)
    return PointGroup(name="I", order=60)
