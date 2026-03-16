"""Adapter for computing dart-map automorphisms via pynauty (nauty).

This module encodes a :class:`~polygraph.structures.dart_map.DartMap` as an
undirected vertex-colored graph that pynauty can analyse, then extracts dart
permutations from pynauty's raw output.

Encoding overview
-----------------
Given *n* darts we build a graph with ``5*n//2`` vertices in three color
classes:

* **Dart vertices** ``0..n-1`` (color 0) — one per dart.
* **Phi-arc gadgets** ``n..2n-1`` (color 1) — gadget ``n+d`` sits between
  dart *d* and ``phi(d)``, encoding the directed face-traversal step.
* **Alpha-arc gadgets** ``2n..2n+n//2-1`` (color 2) — one per undirected
  edge, connecting the two darts of that edge.

An automorphism of this colored graph that permutes dart vertices induces
either ``pi ∘ phi = phi ∘ pi`` (orientation-preserving) or
``pi ∘ phi = phi⁻¹ ∘ pi`` (orientation-reversing), combined with
``pi ∘ alpha = alpha ∘ pi``.  Together these are exactly the full polyhedral
automorphisms, including reflections.

Using phi-arcs (rather than sigma-arcs) is critical: phi-cycles correspond to
faces, so gadgets respect face sizes.  An automorphism can only map a face of
size *k* to another face of size *k*, preventing spurious cross-face-size
symmetries that sigma-arc gadgets would admit.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation

_IMPORT_ERROR_MSG = (
    "pynauty is required for automorphism computation. "
    "Install it with: pip install 'polygraph[symmetry]'"
)


def _build_adjacency(
    n: int,
    phi: list[int],
    phi_inv: list[int],
    alpha: list[int],
    dart_to_alpha_gadget: list[int],
) -> dict[int, list[int]]:
    """Build the adjacency dict for the phi-arc + alpha-arc gadget graph."""
    total = n + n + n // 2
    adj: dict[int, list[int]] = {v: [] for v in range(total)}

    for d in range(n):
        phi_gadget_out = n + d           # gadget for d → phi[d]
        phi_gadget_in = n + phi_inv[d]   # gadget whose target is d
        alpha_gadget = dart_to_alpha_gadget[d]
        adj[d].extend([phi_gadget_out, phi_gadget_in, alpha_gadget])
        adj[phi_gadget_out].extend([d, phi[d]])

    for d in range(n):
        if d < alpha[d]:
            g = dart_to_alpha_gadget[d]
            adj[g].extend([d, alpha[d]])

    return adj


def dartmap_to_pynauty_graph(dm: DartMap) -> object:
    """Build a pynauty Graph encoding the combinatorial structure of *dm*.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    pynauty.Graph
        Vertex-colored undirected graph whose automorphisms correspond to
        dart-map automorphisms.

    Raises
    ------
    ImportError
        If pynauty is not installed.
    """
    try:
        import pynauty  # type: ignore[import-untyped]
    except ImportError as exc:
        raise ImportError(_IMPORT_ERROR_MSG) from exc

    n = dm.num_darts
    alpha = dm.alpha
    phi = [dm.phi(d) for d in range(n)]

    phi_inv: list[int] = [0] * n
    for d, p in enumerate(phi):
        phi_inv[p] = d

    dart_to_alpha_gadget: list[int] = [0] * n
    edge_j = 0
    for d in range(n):
        if d < alpha[d]:
            gadget = 2 * n + edge_j
            dart_to_alpha_gadget[d] = gadget
            dart_to_alpha_gadget[alpha[d]] = gadget
            edge_j += 1

    adj = _build_adjacency(n, phi, phi_inv, alpha, dart_to_alpha_gadget)
    total = n + n + n // 2
    coloring = [
        set(range(n)),
        set(range(n, 2 * n)),
        set(range(2 * n, total)),
    ]

    return pynauty.Graph(
        total,
        directed=False,
        adjacency_dict=adj,
        vertex_coloring=coloring,
    )


def dartmap_automorphisms(
    dm: DartMap,
) -> tuple[list[Permutation], int]:
    """Return the automorphism generators and group order for *dm*.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    generators : list[Permutation]
        Generating set for ``Aut(dm)`` as permutations on dart indices
        ``0..dm.num_darts-1``.
    group_order : int
        ``|Aut(dm)|``.

    Raises
    ------
    ImportError
        If pynauty is not installed.
    """
    try:
        import pynauty  # type: ignore[import-untyped]
    except ImportError as exc:
        raise ImportError(_IMPORT_ERROR_MSG) from exc

    g = dartmap_to_pynauty_graph(dm)
    raw_gens, grpsize1, grpsize2, _orbits, _numorbits = pynauty.autgrp(g)

    n = dm.num_darts
    generators = [Permutation.from_sequence(gen[:n]) for gen in raw_gens]
    group_order = int(round(grpsize1 * (10**grpsize2)))
    return generators, group_order
