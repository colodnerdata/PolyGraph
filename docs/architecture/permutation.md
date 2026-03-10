# Permutation Architecture

## Purpose

`Permutation` is a small, immutable utility type in the `structures` layer for
finite permutations on integer indices.

It provides reusable algebra over index sets without coupling to geometry,
graphs, or polyhedral data structures.

## Scope in the Layered Design

In PolyGraph's architecture (`structures -> algorithms -> geometry ->
visualization`), this module belongs in `structures` because it defines
fundamental combinatorial machinery used by higher-level topology code.

`Permutation` is generic. It can represent:

- vertex rotation maps
- edge pairings
- face successor maps
- symmetry actions on discrete elements

## Representation and Invariants

`Permutation` stores one-line notation as:

- `mapping[i] = image of i`

Validation in `__post_init__` requires:

- `mapping` is a bijection on `0..n-1`, where `n = len(mapping)`

So every index has exactly one image and one preimage.

## Operation Semantics

### Inverse

- `inverse()` returns `p_inv` with `p_inv[p[i]] == i`.

### Composition

- `compose(other)` applies `other` first, then `self`.
- Equivalent index rule:
  - `self.compose(other)[i] == self[other[i]]`
- Both permutations must have the same size.

### Orbit Traversal

- `orbit(start)` yields one full cycle beginning at `start`.
- It stops immediately before `start` would repeat.

### Disjoint Cycles

- `cycles()` decomposes the permutation into disjoint cycles that cover the
  full domain.
- Output order is traversal-based and deterministic for a given mapping:
  scan starts at index `0`, then the next unseen index.
- Cycle presentation is not canonical under relabeling.

## Relationship to `dart_map.py`

`DartMap` is defined by permutations (`sigma`, `alpha`, derived `phi`) over
dart indices. The `Permutation` abstraction models the same mathematical
objects but as a standalone, reusable type.

Current `dart_map.py` operations are implemented directly over `list[int]` for
simplicity, but the semantics match `Permutation` operations exactly:

- inversion
- composition order
- orbit/cycle decomposition

This alignment keeps the topology model consistent and leaves room for future
refactoring toward stronger type-level permutation APIs.

## Minimal Examples

```python
from polygraph.structures.permutation import Permutation

p = Permutation.from_sequence([2, 0, 1])  # 0->2, 1->0, 2->1
q = Permutation.identity(3)

assert p.inverse()[2] == 0
assert p.compose(q).mapping == p.mapping
assert list(p.orbit(0)) == [0, 2, 1]
assert p.cycles() == [[0, 2, 1]]
```
