"""Tests for Johnson solid family generators."""

from __future__ import annotations

import pytest

from polygraph.generators.johnson import (
    dipyramid,
    gyroelongated_square_bipyramid,
    pyramid,
    snub_disphenoid,
    triaugmented_triangular_prism,
)
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
    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())

    assert (num_vertices, num_edges, num_faces) == expected_counts
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
    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())

    assert (num_vertices, num_edges, num_faces) == expected_counts
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


def test_snub_disphenoid_counts() -> None:
    dm = snub_disphenoid()

    assert isinstance(dm, DartMap)
    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())

    assert (num_vertices, num_edges, num_faces) == (8, 18, 12)
    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0


def test_triaugmented_triangular_prism_counts() -> None:
    dm = triaugmented_triangular_prism()

    assert isinstance(dm, DartMap)
    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())

    assert (num_vertices, num_edges, num_faces) == (9, 21, 14)
    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0


def test_gyroelongated_square_bipyramid_counts() -> None:
    dm = gyroelongated_square_bipyramid()

    assert isinstance(dm, DartMap)
    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())

    assert (num_vertices, num_edges, num_faces) == (10, 24, 16)
    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0
