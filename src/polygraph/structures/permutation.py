"""Finite permutations on integer index sets.

This module defines :class:`Permutation`, an immutable permutation type on the
index set ``0..n-1`` stored in one-line notation (``mapping[i]`` is the image
of ``i``).

Composition convention
---------------------
``self.compose(other)`` means "apply ``other`` first, then ``self``", so
``self.compose(other)[i] == self[other[i]]``.

The type is small and independent from geometry. It is intended as a reusable
building block for topology-layer code such as dart-map permutations.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Permutation:
    """Finite permutation of integer indices.

    Parameters
    ----------
    mapping : tuple[int, ...]
        Permutation in one-line notation, where ``mapping[i]`` is the image of
        index ``i``.

    Notes
    -----
    Valid instances represent a bijection on ``0..n-1`` where
    ``n = len(mapping)``.
    """

    mapping: tuple[int, ...]

    def __post_init__(self) -> None:
        """Validate that ``mapping`` is a permutation of ``0..n-1``.

        Raises
        ------
        ValueError
            If ``mapping`` is not a bijection on ``0..n-1``.
        """
        n = len(self.mapping)
        if set(self.mapping) != set(range(n)):
            raise ValueError("mapping must be a permutation of 0..n-1")

    def __len__(self) -> int:
        """Return the size of the permutation domain.

        Returns
        -------
        int
            Number of elements permuted.
        """
        return len(self.mapping)

    def __getitem__(self, i: int) -> int:
        """Return the image of an index under the permutation.

        Parameters
        ----------
        i : int
            Index to evaluate. Standard tuple indexing semantics apply.

        Returns
        -------
        int
            Image of ``i`` under the permutation.

        Raises
        ------
        IndexError
            If ``i`` is outside valid tuple index bounds.
        """
        return self.mapping[i]

    def inverse(self) -> Permutation:
        """Return the inverse permutation.

        Returns
        -------
        Permutation
            Permutation ``inv`` such that ``inv[self[i]] == i`` for all valid
            indices ``i``.
        """
        inv = [0] * len(self.mapping)
        for i, j in enumerate(self.mapping):
            inv[j] = i
        return Permutation(tuple(inv))

    def compose(self, other: Permutation) -> Permutation:
        """Compose two permutations of the same size.

        Parameters
        ----------
        other : Permutation
            Permutation to apply first.

        Returns
        -------
        Permutation
            Composition ``self.compose(other)`` where
            ``result[i] == self[other[i]]``.

        Raises
        ------
        ValueError
            If ``self`` and ``other`` do not have the same size.

        Notes
        -----
        Composition order follows function composition: ``self`` after
        ``other``.
        """
        if len(self) != len(other):
            raise ValueError("permutations must have the same size")
        return Permutation(tuple(self[other[i]] for i in range(len(self))))

    def orbit(self, start: int) -> Iterator[int]:
        """Yield one full orbit under repeated application.

        Parameters
        ----------
        start : int
            Starting index for orbit traversal.

        Yields
        ------
        int
            Successive elements in the orbit beginning at ``start`` and ending
            immediately before repeating ``start``.

        Raises
        ------
        IndexError
            If ``start`` is invalid when the permutation is applied.

        Notes
        -----
        This generator emits each orbit element exactly once in traversal order.
        """
        i = start
        while True:
            yield i
            i = self[i]
            if i == start:
                break

    def cycles(self) -> list[list[int]]:
        """Return disjoint cycle decomposition.

        Returns
        -------
        list[list[int]]
            Disjoint cycles covering all indices in ``0..n-1``.

        Notes
        -----
        Cycle order is traversal-based: cycles are discovered by scanning
        indices from ``0`` upward, and each cycle follows repeated application
        of the permutation. The output is not canonical under relabeling.
        """
        seen = [False] * len(self.mapping)
        out: list[list[int]] = []

        for i in range(len(self.mapping)):
            if seen[i]:
                continue
            cycle = []
            j = i
            while not seen[j]:
                seen[j] = True
                cycle.append(j)
                j = self[j]
            out.append(cycle)

        return out

    @classmethod
    def identity(cls, n: int) -> Permutation:
        """Construct the identity permutation.

        Parameters
        ----------
        n : int
            Domain size.

        Returns
        -------
        Permutation
            Identity permutation on ``0..n-1``.
        """
        return cls(tuple(range(n)))

    @classmethod
    def from_sequence(cls, mapping: Sequence[int]) -> Permutation:
        """Construct a permutation from a sequence.

        Parameters
        ----------
        mapping : Sequence[int]
            Candidate permutation in one-line notation.

        Returns
        -------
        Permutation
            New immutable permutation instance.

        Raises
        ------
        ValueError
            If ``mapping`` is not a permutation of ``0..n-1``.
        """
        return cls(tuple(mapping))
