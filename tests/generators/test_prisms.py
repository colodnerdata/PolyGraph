"""Tests for prism and antiprism generators."""

from __future__ import annotations

import pytest

from polygraph.generators.prisms import antiprism, prism


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_prism_returns_sphere_map(n: int) -> None:
    """Prism generator should produce closed genus-zero maps for various n.

    For an n-gonal prism:
    - vertices = 2n
    - edges = 3n  (n top + n bottom + n lateral)
    - faces = n + 2  (n side quads + 2 caps)
    """
    dm = prism(n)

    assert len(dm.vertex_orbits()) == 2 * n
    assert dm.num_edges == 3 * n
    assert len(dm.face_orbits()) == n + 2
    assert len(dm.connected_components()) == 1
    assert dm.genus() == 0


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_antiprism_returns_sphere_map(n: int) -> None:
    """Antiprism generator should produce closed genus-zero maps for various n.

    For an n-gonal antiprism:
    - vertices = 2n
    - edges = 4n  (n top + n bottom + 2n lateral)
    - faces = 2n + 2  (2n triangles + 2 caps)
    """
    dm = antiprism(n)

    assert len(dm.vertex_orbits()) == 2 * n
    assert dm.num_edges == 4 * n
    assert len(dm.face_orbits()) == 2 * n + 2
    assert len(dm.connected_components()) == 1
    assert dm.genus() == 0


@pytest.mark.parametrize("factory", [prism, antiprism])
def test_prism_family_generators_require_at_least_three_sides(factory) -> None:
    """Prism-family generators should reject degenerate polygon counts."""
    with pytest.raises(ValueError, match=r"requires n >= 3"):
        factory(2)
