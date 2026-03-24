"""Shared pytest fixtures for DartMap test suites."""

from __future__ import annotations

import pytest

from polygraph.generators.platonic import (
    cube,
    dodecahedron,
    icosahedron,
    octahedron,
    tetrahedron,
)


@pytest.fixture(scope="session")
def tetrahedron_dm():
    return tetrahedron()


@pytest.fixture(scope="session")
def cube_dm():
    return cube()


@pytest.fixture(scope="session")
def octahedron_dm():
    return octahedron()


@pytest.fixture(scope="session")
def dodecahedron_dm():
    return dodecahedron()


@pytest.fixture(scope="session")
def icosahedron_dm():
    return icosahedron()


@pytest.fixture(
    scope="session",
    params=[
        "tetrahedron_dm",
        "cube_dm",
        "octahedron_dm",
        "dodecahedron_dm",
        "icosahedron_dm",
    ],
)
def platonic_dm(request):
    return request.getfixturevalue(request.param)
