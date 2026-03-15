"""pynauty-based adapter for DartMap automorphism computation.

Encodes a DartMap as an **undirected**, vertex-colored graph for pynauty's
canonical labeling engine (nauty).  The encoding uses:

- **Dart vertices** (indices ``0..n-1``, color 0): one per dart.
- **Phi-arc vertices** (indices ``n..2n-1``, color 1): one per dart ``d``,
  connecting dart ``d`` to dart ``phi[d] = sigma⁻¹(alpha[d])``
  via an undirected edge ``d — n+d — phi[d]``.
- **Alpha-arc vertices** (indices ``2n..2n+n//2-1``, color 2): one per
  *undirected* edge ``{d, alpha[d]}``, connecting dart ``d`` to dart
  ``alpha[d]``.  Only the canonical representative ``d < alpha[d]`` is
  assigned a gadget to avoid duplicates.

**Why phi-arcs instead of sigma-arcs:**
Undirected sigma-arc gadgets admit spurious automorphisms satisfying
``pi ∘ sigma = sigma⁻¹ ∘ pi``.  Such maps send phi-faces to
``sigma ∘ alpha`` orbits — which are not faces — so they are NOT valid
polyhedral symmetries.  Undirected phi-arc gadgets instead admit maps
satisfying ``pi ∘ phi = phi ∘ pi`` (orientation-preserving) or
``pi ∘ phi = phi⁻¹ ∘ pi`` (orientation-reversing, faces map to faces in
the opposite traversal order).  Together with ``pi ∘ alpha = alpha ∘ pi``
(enforced by the alpha gadgets), these are exactly the full polyhedral
automorphism group including reflections.

For connected orientable surfaces (all valid polyhedra), a direction
reversal at any phi-arc propagates to all others via the ``phi``/``alpha``
alternating path, so globally mixed orientation is impossible.

Only ``pynauty`` is used; no BLISS binding is required.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation


def _build_pynauty_graph(dm: DartMap):
    """Build an undirected vertex-colored pynauty Graph from a DartMap.

    Parameters
    ----------
    dm : DartMap
        Input combinatorial map.

    Returns
    -------
    pynauty.Graph
        Undirected graph with ``5 * dm.num_darts // 2`` vertices
        encoding ``phi`` (one gadget per dart arc) and ``alpha``
        (one gadget per undirected edge).

    Raises
    ------
    ImportError
        If ``pynauty`` is not installed.
    """
    try:
        import pynauty
    except ImportError as exc:
        raise ImportError(
            "pynauty is required for symmetry computations. "
            "Install with: pip install 'polygraph[symmetry]'"
        ) from exc

    n = dm.num_darts
    alpha_start = 2 * n
    num_edges = n // 2

    adj: dict[int, list[int]] = {}

    # Phi-arc gadgets: vertex n+d sits between dart d and dart phi(d).
    for d in range(n):
        phi_d = dm._sigma_inv[dm.alpha[d]]
        adj[n + d] = [d, phi_d]

    # Alpha-arc gadgets: one per canonical undirected edge (d < alpha[d]).
    edge_count = 0
    for d in range(n):
        if d < dm.alpha[d]:
            g = alpha_start + edge_count
            adj[g] = [d, dm.alpha[d]]
            edge_count += 1

    coloring = [
        set(range(n)),                                     # dart vertices
        set(range(n, 2 * n)),                              # phi-arc vertices
        set(range(alpha_start, alpha_start + num_edges)),  # alpha-arc vertices
    ]

    return pynauty.Graph(
        number_of_vertices=2 * n + num_edges,
        directed=False,
        adjacency_dict=adj,
        vertex_coloring=coloring,
    )


def dartmap_automorphisms(
    dm: DartMap,
) -> tuple[list[Permutation], int]:
    """Compute automorphism generators of a DartMap using pynauty.

    Computes the full automorphism group of the combinatorial map,
    including both orientation-preserving automorphisms (commuting with
    ``phi``) and orientation-reversing automorphisms (reversing ``phi``).
    For the Platonic solids this yields the full symmetry group orders
    (T_d = 24, O_h = 48, I_h = 120).

    Parameters
    ----------
    dm : DartMap
        Input combinatorial map.

    Returns
    -------
    generators : list[Permutation]
        Generating set for ``Aut(dm)`` as permutations on dart indices
        ``0..dm.num_darts-1``.
    group_order : int
        Order of the automorphism group ``|Aut(dm)|``.

    Raises
    ------
    ImportError
        If ``pynauty`` is not installed.

    Notes
    -----
    The group order is derived from pynauty's ``grpsize1 * 10^grpsize2``
    representation.  For the polyhedral sizes this library targets the
    result is exact.
    """
    import pynauty

    g = _build_pynauty_graph(dm)
    gens_raw, grpsize1, grpsize2, _, _ = pynauty.autgrp(g)

    n = dm.num_darts
    # Project each (5n/2)-element permutation down to dart indices.
    generators = [
        Permutation.from_sequence([gen[d] for d in range(n)])
        for gen in gens_raw
    ]
    group_order = int(round(grpsize1)) * (10 ** int(grpsize2))
    return generators, group_order
