"""Tests for prism and antiprism generators."""

from __future__ import annotations

import pytest

from polygraph.generators.prisms import antiprism, prism


@pytest.mark.parametrize(
    ("factory", "n", "num_edges", "num_faces"),
    [
        (prism, 5, 15, 7),
        (antiprism, 5, 20, 12),
    ],
)
def test_prism_family_generators_return_sphere_maps(
    factory, n: int, num_edges: int, num_faces: int
) -> None:
    """Prism-family generators should produce closed genus-zero maps."""
    dm = factory(n)

    assert len(dm.vertex_orbits()) == 2 * n
    assert dm.num_edges == num_edges
    assert len(dm.face_orbits()) == num_faces
    assert dm.genus() == 0


@pytest.mark.parametrize("factory", [prism, antiprism])
def test_prism_family_generators_require_at_least_three_sides(factory) -> None:
    """Prism-family generators should reject degenerate polygon counts."""
    with pytest.raises(ValueError, match=r"requires n >= 3"):
        factory(2)
