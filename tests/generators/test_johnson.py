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
def test_dipyramid_counts(
    n: int, expected_counts: tuple[int, int, int]
) -> None:
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


# ---------------------------------------------------------------------------
# Johnson deltahedra: snub disphenoid (J84), triaugmented triangular prism
# (J51), and gyroelongated square bipyramid (J17)
# ---------------------------------------------------------------------------


def _assert_closed_map(dm: DartMap) -> None:
    """Assert that alpha is a fixed-point-free involution (closed map)."""
    for d in range(dm.num_darts):
        twin = dm.alpha[d]
        assert twin != d, f"dart {d} is its own twin"
        assert dm.alpha[twin] == d, f"alpha is not an involution at dart {d}"


@pytest.mark.parametrize(
    ("generator", "expected_counts"),
    [
        (snub_disphenoid, (8, 18, 12)),
        (triaugmented_triangular_prism, (9, 21, 14)),
        (gyroelongated_square_bipyramid, (10, 24, 16)),
    ],
)
def test_johnson_deltahedra_counts(
    generator, expected_counts: tuple[int, int, int]
) -> None:
    """New Johnson generators produce maps with the expected (V, E, F)."""
    dm = generator()

    assert isinstance(dm, DartMap)
    num_vertices = len(dm.vertex_orbits())
    num_edges = dm.num_edges
    num_faces = len(dm.face_orbits())

    assert (num_vertices, num_edges, num_faces) == expected_counts


@pytest.mark.parametrize(
    "generator",
    [
        snub_disphenoid,
        triaugmented_triangular_prism,
        gyroelongated_square_bipyramid,
    ],
)
def test_johnson_deltahedra_euler_characteristic(generator) -> None:
    """New Johnson generators produce sphere-topology maps (χ=2, genus=0)."""
    dm = generator()

    assert dm.euler_characteristic() == 2
    assert dm.genus() == 0


@pytest.mark.parametrize(
    "generator",
    [
        snub_disphenoid,
        triaugmented_triangular_prism,
        gyroelongated_square_bipyramid,
    ],
)
def test_johnson_deltahedra_closed_map(generator) -> None:
    """New Johnson generators produce closed maps (no edge-pairing errors)."""
    dm = generator()

    _assert_closed_map(dm)
    assert len(dm.connected_components()) == 1


@pytest.mark.parametrize(
    "generator",
    [
        snub_disphenoid,
        triaugmented_triangular_prism,
        gyroelongated_square_bipyramid,
    ],
)
def test_johnson_deltahedra_all_triangular_faces(generator) -> None:
    """All faces of the Johnson deltahedra generators are triangles."""
    dm = generator()

    face_sizes = {len(face) for face in dm.face_orbits()}
    assert face_sizes == {3}
