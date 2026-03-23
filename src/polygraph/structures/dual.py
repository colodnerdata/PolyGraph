"""Dual-map construction and map-isomorphism checks.

The dual of a combinatorial map ``(sigma, alpha)`` keeps the same dart set and
edge involution ``alpha``, while replacing ``sigma`` with
``phi = sigma^{-1} o alpha``.
"""

from __future__ import annotations

from collections.abc import Sequence

from polygraph.structures.dart_map import DartMap


def dual_of(dm: DartMap) -> DartMap:
    """Return the dual map of ``dm``.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    DartMap
        Dual map ``(phi, alpha)`` on the same darts.
    """
    phi_list = [dm.phi(d) for d in range(dm.num_darts)]
    return DartMap(sigma=phi_list, alpha=list(dm.alpha))


def _invert(perm: Sequence[int], n: int) -> list[int] | None:
    """Return inverse permutation, or ``None`` if ``perm`` is invalid."""
    if len(perm) != n:
        return None

    inv = [-1] * n
    for i, j in enumerate(perm):
        if not 0 <= j < n:
            return None
        if inv[j] != -1:
            return None
        inv[j] = i
    return inv


def is_isomorphism(dm1: DartMap, dm2: DartMap, pi: Sequence[int]) -> bool:
    """Check whether ``pi`` is a dart relabeling isomorphism ``dm1 -> dm2``.

    The map must conjugate both structure permutations:

    - ``pi * dm1.sigma * pi^{-1} == dm2.sigma``
    - ``pi * dm1.alpha * pi^{-1} == dm2.alpha``

    Invalid ``pi`` values (wrong length or non-permutations) return ``False``.
    """
    if dm1.num_darts != dm2.num_darts:
        return False

    n = dm1.num_darts
    pi_inv = _invert(pi, n)
    if pi_inv is None:
        return False

    for d in range(n):
        if pi[dm1.sigma[pi_inv[d]]] != dm2.sigma[d]:
            return False
        if pi[dm1.alpha[pi_inv[d]]] != dm2.alpha[d]:
            return False
    return True


__all__ = ["dual_of", "is_isomorphism"]
