"""Tests for Johnson solid family generators."""

from __future__ import annotations

import pytest

from polygraph.generators.johnson import dipyramid, pyramid
from polygraph.structures.dart_map import DartMap


@pytest.mark.parametrize(
    ("n", "expected_counts"),
    [
        (3, (4, 6, 4)),
        (4, (5, 8, 5)),
        (5, (6, 10, 6)),
    ],
)
def test_pyramid_counts(n: int, expected_counts: tuple[int, int, int]) -> None:
    dm = pyramid(n)

    assert isinstance(dm, DartMap)
    assert (len(dm.vertex_orbits()), dm.num_edges, len(dm.face_orbits())) == expected_counts
    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0


@pytest.mark.parametrize(
    ("n", "expected_counts"),
    [
        (3, (5, 9, 6)),
        (4, (6, 12, 8)),
        (5, (7, 15, 10)),
    ],
)
def test_dipyramid_counts(n: int, expected_counts: tuple[int, int, int]) -> None:
    dm = dipyramid(n)

    assert isinstance(dm, DartMap)
    assert (len(dm.vertex_orbits()), dm.num_edges, len(dm.face_orbits())) == expected_counts
    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0


@pytest.mark.parametrize("n", [0, 1, 2])
def test_pyramid_requires_at_least_three_sides(n: int) -> None:
    with pytest.raises(ValueError, match=r"requires n >= 3"):
        pyramid(n)


@pytest.mark.parametrize("n", [0, 1, 2])
def test_dipyramid_requires_at_least_three_sides(n: int) -> None:
    with pytest.raises(ValueError, match=r"requires n >= 3"):
        dipyramid(n)
