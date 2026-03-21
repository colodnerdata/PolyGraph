"""Tests for Schlaefli and vertex-configuration notation helpers."""

from __future__ import annotations

import pytest

from polygraph.generators.johnson import pyramid
from polygraph.generators.notation import schlafli_symbol, vertex_configuration
from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)
from polygraph.generators.prisms import antiprism, prism

# ---------------------------------------------------------------------------
# schlafli_symbol
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("generator", "expected"),
    [
        (tetrahedron, "{3, 3}"),
        (cube, "{4, 3}"),
        (octahedron, "{3, 4}"),
        (dodecahedron, "{5, 3}"),
        (icosahedron, "{3, 5}"),
    ],
)
def test_schlafli_symbol_returns_expected_pair(
    generator,
    expected: str,
) -> None:
    assert schlafli_symbol(generator()) == expected


def test_schlafli_symbol_rejects_non_regular_polyhedron() -> None:
    with pytest.raises(ValueError, match=r"Not a regular polyhedron"):
        schlafli_symbol(prism(3))


# ---------------------------------------------------------------------------
# vertex_configuration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("dm", "expected"),
    [
        (tetrahedron(), "3.3.3"),
        (cube(), "4.4.4"),
        (octahedron(), "3.3.3.3"),
        (dodecahedron(), "5.5.5"),
        (icosahedron(), "3.3.3.3.3"),
        (prism(3), "3.4.4"),
        (prism(5), "4.4.5"),
        (prism(6), "4.4.6"),
        (antiprism(3), "3.3.3.3"),
        (antiprism(4), "3.3.3.4"),
        (antiprism(5), "3.3.3.5"),
    ],
)
def test_vertex_configuration_returns_expected_string(
    dm,
    expected: str,
) -> None:
    assert vertex_configuration(dm) == expected


def test_vertex_configuration_rejects_non_vertex_transitive_polyhedron(
) -> None:
    with pytest.raises(ValueError, match=r"Not all vertices share the same"):
        vertex_configuration(pyramid(4))


@pytest.mark.parametrize("n", [4, 5, 6])
def test_vertex_configuration_rejects_all_pyramids(n: int) -> None:
    with pytest.raises(ValueError):
        vertex_configuration(pyramid(n))


# ---------------------------------------------------------------------------
# Consistency: schlafli p matches every face size in vertex_configuration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [tetrahedron, cube, octahedron, dodecahedron, icosahedron],
)
def test_schlafli_p_matches_vertex_config(factory) -> None:
    """For regular solids all face sizes in vc equal p."""
    dm = factory()
    symbol = schlafli_symbol(dm)
    p = int(symbol.removeprefix("{").removesuffix("}").split(",")[0])
    assert all(int(s) == p for s in vertex_configuration(dm).split("."))