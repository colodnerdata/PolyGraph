"""Validation utilities for triangulated dart maps."""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap


def validate_triangulation(original: DartMap, subdivided: DartMap) -> None:
    """Validate that a subdivided map is a correct barycentric subdivision.

    Parameters
    ----------
    original : DartMap
        The original dart map before subdivision.
    subdivided : DartMap
        The subdivided dart map to validate.

    Raises
    ------
    ValueError
        If any validation check fails.
    """
    n = original.num_darts
    e_orig = original.num_edges
    v_orig = len(original.vertex_orbits())
    f_orig = len(original.face_orbits())

    # Check dart count: 6n
    expected_darts = 6 * n
    if subdivided.num_darts != expected_darts:
        raise ValueError(
            f"Expected {expected_darts} darts in subdivision, "
            f"got {subdivided.num_darts}."
        )

    # Check all faces are triangles (phi orbits of length 3)
    for face in subdivided.face_orbits():
        if len(face) != 3:
            raise ValueError(
                f"Subdivision face has {len(face)} darts; expected 3."
            )

    # Check face count: 4E = 2n
    face_count = len(subdivided.face_orbits())
    expected_faces = 2 * n
    if face_count != expected_faces:
        raise ValueError(
            f"Expected {expected_faces} faces, got {face_count}."
        )

    # Check vertex count: V + E + F
    vertex_count = len(subdivided.vertex_orbits())
    expected_vertices = v_orig + e_orig + f_orig
    if vertex_count != expected_vertices:
        raise ValueError(
            f"Expected {expected_vertices} vertices (V+E+F = "
            f"{v_orig}+{e_orig}+{f_orig}), got {vertex_count}."
        )

    # Check edge count: 6E = 3n
    edge_count = subdivided.num_edges
    expected_edges = 3 * n
    if edge_count != expected_edges:
        raise ValueError(
            f"Expected {expected_edges} edges, got {edge_count}."
        )

    # Check Euler characteristic preserved
    chi_orig = original.euler_characteristic()
    chi_sub = subdivided.euler_characteristic()
    if chi_sub != chi_orig:
        raise ValueError(
            f"Euler characteristic changed: {chi_orig} → {chi_sub}."
        )

    # Check genus preserved
    g_orig = original.genus()
    g_sub = subdivided.genus()
    if g_sub != g_orig:
        raise ValueError(
            f"Genus changed: {g_orig} → {g_sub}."
        )
