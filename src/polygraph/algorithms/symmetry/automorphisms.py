"""Automorphism group computation for combinatorial maps.

The primary computation path uses ``pynauty`` (nauty backend) via
:func:`polygraph.interop.bliss_adapter.dartmap_automorphisms`.

Public API
----------
compute_automorphism_generators
    Return a generating set for ``Aut(dm)``.
automorphism_group_order
    Compute ``|Aut(dm)|`` from a generator set via BFS group closure.
is_orientation_preserving
    Check whether an automorphism commutes with ``sigma``.
"""

from __future__ import annotations

from polygraph.interop.bliss_adapter import dartmap_automorphisms
from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation

# Maximum group order for which BFS group closure is attempted.
_BFS_ORDER_LIMIT = 100_000


def compute_automorphism_generators(dm: DartMap) -> list[Permutation]:
    """Return a generating set for the automorphism group of a DartMap.

    Parameters
    ----------
    dm : DartMap
        Input combinatorial map.

    Returns
    -------
    list[Permutation]
        Generators of ``Aut(dm)`` as permutations on dart indices.
        The list may be empty when the map has no non-trivial symmetry.

    Raises
    ------
    ImportError
        If ``pynauty`` is not installed.
    """
    generators, _ = dartmap_automorphisms(dm)
    return generators


def automorphism_group_order(
    generators: list[Permutation], num_darts: int
) -> int:
    """Compute the automorphism group order from a generator set.

    Enumerates all group elements reachable from the identity by
    repeated application of the generators and their inverses.  This
    BFS closure approach is exact for groups up to
    ``_BFS_ORDER_LIMIT`` elements; larger groups raise ``ValueError``.

    Parameters
    ----------
    generators : list[Permutation]
        Generating set for the group, as returned by
        :func:`compute_automorphism_generators`.
    num_darts : int
        Domain size (number of darts in the map).

    Returns
    -------
    int
        ``|Aut(dm)|``, including the identity.

    Raises
    ------
    ValueError
        If group enumeration exceeds ``_BFS_ORDER_LIMIT`` elements.

    Notes
    -----
    For maps with large symmetry groups (e.g., Goldberg polyhedra),
    prefer the ``group_order`` value returned directly by
    :func:`compute_automorphism_generators` via
    :func:`~polygraph.interop.bliss_adapter.dartmap_automorphisms`.
    """
    if not generators:
        return 1

    identity = Permutation.identity(num_darts)
    # Include inverses so a single BFS pass suffices.
    all_gens = generators + [g.inverse() for g in generators]

    seen: set[Permutation] = {identity}
    queue: list[Permutation] = [identity]

    while queue:
        current = queue.pop()
        for gen in all_gens:
            product = gen.compose(current)
            if product not in seen:
                seen.add(product)
                if len(seen) > _BFS_ORDER_LIMIT:
                    raise ValueError(
                        f"Group order exceeds BFS limit of "
                        f"{_BFS_ORDER_LIMIT}.  Use the group_order "
                        "value from dartmap_automorphisms() instead."
                    )
                queue.append(product)

    return len(seen)


def is_orientation_preserving(pi: Permutation, dm: DartMap) -> bool:
    """Return whether an automorphism preserves face orientation.

    An automorphism ``pi`` is orientation-preserving when it commutes
    with the vertex-rotation permutation ``sigma``:

        ``pi(sigma(d)) == sigma(pi(d))``  for all darts ``d``.

    All automorphisms returned by :func:`compute_automorphism_generators`
    satisfy this condition by construction.  This function is provided
    for validation and for detecting orientation-reversing symmetries if
    they are supplied by the caller.

    Parameters
    ----------
    pi : Permutation
        A permutation on dart indices ``0..dm.num_darts-1``.
    dm : DartMap
        The map ``pi`` is claimed to act on.

    Returns
    -------
    bool
        ``True`` iff ``pi`` commutes with ``sigma``.

    Raises
    ------
    ValueError
        If the domain size of ``pi`` does not match ``dm.num_darts``.
    """
    if len(pi) != dm.num_darts:
        raise ValueError(
            f"Permutation domain size {len(pi)} does not match "
            f"dm.num_darts = {dm.num_darts}."
        )
    sigma = dm.sigma
    for d in range(dm.num_darts):
        if pi[sigma[d]] != sigma[pi[d]]:
            return False
    return True
