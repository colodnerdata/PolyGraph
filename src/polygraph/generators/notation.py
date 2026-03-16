"""Schläfli symbols and vertex-configuration strings for polyhedral dart maps.

Both functions are pure combinatorial queries on the dart map — no geometry
is required.  They rely only on the traversal layer.

**Schläfli symbol** ``{p, q}``
    Defined only for *regular* polyhedra where every face is a *p*-gon and
    every vertex has degree *q*.  Raises :exc:`ValueError` for non-regular
    inputs.

**Vertex-configuration string**
    Generalises the Schläfli symbol to *semi-regular* (Archimedean, prismatic)
    polyhedra.  Records the ordered sequence of face sizes around each vertex
    in sigma-rotation order, e.g. ``"3.4.4"`` for a triangular prism.
    Raises :exc:`ValueError` if vertex configurations differ (not
    vertex-transitive).
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.traversal import (
    all_face_orbits,
    all_vertex_orbits,
    face_darts,
    vertex_darts,
)

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _min_rotation(seq: tuple[int, ...]) -> tuple[int, ...]:
    """Return the lexicographically smallest cyclic rotation of *seq*."""
    if not seq:
        return ()
    n = len(seq)
    return min(seq[i:] + seq[:i] for i in range(n))


def _face_sizes_at_vertex(dm: DartMap, v: int) -> tuple[int, ...]:
    """Return face sizes around one vertex in sigma-rotation order.

    Parameters
    ----------
    dm : DartMap
        Input dart map.
    v : int
        Representative dart of the vertex orbit to inspect.

    Returns
    -------
    face_sizes : tuple[int, ...]
        One integer per incident face — the number of boundary darts on that
        face — in the order they are encountered walking sigma around *v*.

    Notes
    -----
    This is a pure combinatorial query: for each dart in the vertex orbit of
    ``v``, it counts the darts in the corresponding face orbit.

    Implementation details
    ----------------------
    Directly calling :func:`face_darts` for every incident dart would cause
    the underlying ``phi`` permutation to be rebuilt repeatedly.  Instead we
    precompute the size of every face once and map each dart to the size of
    its incident face, then perform cheap lookups while walking around the
    vertex.
    """
    # Precompute face sizes for all darts in the map.
    face_size_by_dart: dict[int, int] = {}
    for f in all_face_orbits(dm):
        darts_on_face = tuple(face_darts(dm, f))
        size = len(darts_on_face)
        for d in darts_on_face:
            face_size_by_dart[d] = size

    # Collect face sizes around the given vertex in sigma-rotation order.
    return tuple(face_size_by_dart[d] for d in vertex_darts(dm, v))

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def schlafli_symbol(dm: DartMap) -> str:
    """Return the Schläfli symbol ``{p, q}`` for a regular polyhedron.

    Parameters
    ----------
    dm : DartMap
        Input dart map.  Must be a regular polyhedron (all faces congruent,
        all vertices equivalent).

    Returns
    -------
    str
        Schlaefli symbol formatted as ``"{p, q}"`` where every face is a
        *p*-gon and every vertex has degree *q*.

    Raises
    ------
    ValueError
        If face sizes or vertex degrees are not uniform.

    Examples
    --------
    >>> from polygraph.generators.platonic import cube
    >>> schlafli_symbol(cube())
    '{4, 3}'
    """
    face_sizes = {
        sum(1 for _ in face_darts(dm, f)) for f in all_face_orbits(dm)
    }
    vertex_degrees = {
        sum(1 for _ in vertex_darts(dm, v)) for v in all_vertex_orbits(dm)
    }
    if len(face_sizes) != 1 or len(vertex_degrees) != 1:
        raise ValueError(
            "Not a regular polyhedron: "
            f"face sizes={sorted(face_sizes)}, "
            f"vertex degrees={sorted(vertex_degrees)}"
        )
    p = face_sizes.pop()
    q = vertex_degrees.pop()
    return f"{{{p}, {q}}}"


def vertex_configuration(dm: DartMap) -> str:
    """Return the canonical vertex-configuration string.

    Parameters
    ----------
    dm : DartMap
        Input dart map.  Must be vertex-transitive (every vertex has the same
        sequence of surrounding face sizes up to cyclic rotation).

    Returns
    -------
    str
        Dot-separated face sizes, canonicalised to the lexicographically
        smallest cyclic rotation, e.g. ``"3.4.4"`` for a triangular prism,
        ``"3.3.3.3.3"`` for the icosahedron.

    Raises
    ------
    ValueError
        If vertex configurations differ between vertices.

    Examples
    --------
    >>> from polygraph.generators.prisms import prism
    >>> vertex_configuration(prism(3))
    '3.4.4'
    """
    vertices = list(all_vertex_orbits(dm))
    configs = [_face_sizes_at_vertex(dm, v) for v in vertices]
    canonical = [_min_rotation(c) for c in configs]

    if len(set(canonical)) != 1:
        raise ValueError(
            "Polyhedron is not vertex-transitive: "
            "vertex configurations differ between vertices"
        )
    return ".".join(str(s) for s in canonical[0])
