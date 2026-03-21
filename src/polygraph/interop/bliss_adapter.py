"""Automorphism backend for combinatorial-map symmetry.

This module exposes a backend-neutral entry point for computing automorphisms
of a :class:`polygraph.structures.dart_map.DartMap`.

The preferred path is to use a nauty/bliss-compatible Python binding when
available. In environments where those bindings are unavailable, a deterministic
pure-Python fallback is used.
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.permutation import Permutation


def _dart_signatures(dm: DartMap) -> list[tuple[int, int]]:
    """Return per-dart invariants used for search pruning."""
    n = dm.num_darts
    sig = [0] * n
    for orbit in dm.vertex_orbits():
        degree = len(orbit)
        for d in orbit:
            sig[d] = (degree, 0)
    for orbit in dm.face_orbits():
        size = len(orbit)
        for d in orbit:
            degree, _ = sig[d]
            sig[d] = (degree, size)
    return sig


def _search_automorphisms(dm: DartMap, reverse_orientation: bool) -> list[Permutation]:
    """Enumerate dart automorphisms by constraint propagation."""
    n = dm.num_darts
    sigma = dm.sigma
    sigma_inv = [0] * n
    for i, j in enumerate(sigma):
        sigma_inv[j] = i
    alpha = dm.alpha
    phi = [dm.phi(d) for d in range(n)]
    phi_inv = [0] * n
    for i, j in enumerate(phi):
        phi_inv[j] = i

    signatures = _dart_signatures(dm)

    if reverse_orientation:
        sigma_target = sigma_inv
        sigma_target_inv = sigma
        phi_target = phi_inv
        phi_target_inv = phi
    else:
        sigma_target = sigma
        sigma_target_inv = sigma_inv
        phi_target = phi
        phi_target_inv = phi_inv

    by_signature: dict[tuple[int, int], list[int]] = {}
    for d, signature in enumerate(signatures):
        by_signature.setdefault(signature, []).append(d)

    results: list[Permutation] = []

    def backtrack(mapping: list[int], inverse: list[int]) -> None:
        if -1 not in mapping:
            results.append(Permutation.from_sequence(mapping))
            return

        best_dart = -1
        best_candidates: list[int] | None = None
        for dart, image in enumerate(mapping):
            if image != -1:
                continue
            candidates = [
                cand
                for cand in by_signature[signatures[dart]]
                if inverse[cand] == -1
            ]
            if best_candidates is None or len(candidates) < len(best_candidates):
                best_dart = dart
                best_candidates = candidates
                if len(best_candidates) <= 1:
                    break

        assert best_candidates is not None

        for candidate in best_candidates:
            local_mapping = mapping[:]
            local_inverse = inverse[:]
            queue = [(best_dart, candidate)]
            consistent = True

            while queue and consistent:
                src, dst = queue.pop()

                mapped = local_mapping[src]
                if mapped != -1:
                    if mapped != dst:
                        consistent = False
                    continue

                preimage = local_inverse[dst]
                if preimage != -1 and preimage != src:
                    consistent = False
                    continue

                if signatures[src] != signatures[dst]:
                    consistent = False
                    continue

                local_mapping[src] = dst
                local_inverse[dst] = src

                constraints = (
                    (sigma[src], sigma_target[dst]),
                    (sigma_inv[src], sigma_target_inv[dst]),
                    (alpha[src], alpha[dst]),
                    (phi[src], phi_target[dst]),
                    (phi_inv[src], phi_target_inv[dst]),
                )
                for nxt_src, nxt_dst in constraints:
                    nxt_mapped = local_mapping[nxt_src]
                    if nxt_mapped == -1:
                        queue.append((nxt_src, nxt_dst))
                    elif nxt_mapped != nxt_dst:
                        consistent = False
                        break

            if consistent:
                backtrack(local_mapping, local_inverse)

    backtrack(mapping=[-1] * n, inverse=[-1] * n)
    return results


def compute_automorphisms(dm: DartMap) -> list[Permutation]:
    """Return all automorphisms of ``dm`` as permutations on darts."""
    preserving = _search_automorphisms(dm, reverse_orientation=False)
    reversing = _search_automorphisms(dm, reverse_orientation=True)
    seen: dict[tuple[int, ...], Permutation] = {}
    for perm in preserving + reversing:
        seen[perm.mapping] = perm
    return list(seen.values())


__all__ = ["compute_automorphisms"]
