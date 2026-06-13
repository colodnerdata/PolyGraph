"""Planar embedding view over a DartMap.

Provides integer-indexed vertex and face queries needed by planar drawing
algorithms.  Vertex IDs ``0..V-1`` and face IDs ``0..F-1`` are assigned in
``vertex_orbits()`` / ``face_orbits()`` scan order — the same convention used
by :mod:`polygraph.interop.networkx_adapter`.

Orientation notes
-----------------
``sigma`` cycles darts **counter-clockwise** around each vertex in the plane
(verified empirically: when interior faces are wound CCW, the sigma orbit at
a vertex visits neighbours in CCW angular order).  Consequently:

- :meth:`ordered_neighbors` returns neighbour IDs in **CCW** order.
- :meth:`face_boundary_vertices` returns vertices in the **phi orbit** order,
  which is **CCW** for interior faces and **CW** for the outer face (when the
  outer face is listed CW in ``from_face_lists``).
"""

from __future__ import annotations

from polygraph.structures.dart_map import DartMap
from polygraph.structures.traversal import (
    face_darts,
    vertex_darts,
)


def _build_dart_to_vertex(dm: DartMap) -> list[int]:
    """Build a dart-index → vertex-ID lookup array.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[int]
        Array of length ``dm.num_darts`` mapping each dart to its vertex ID.
    """
    dart_to_vertex: list[int] = [-1] * dm.num_darts
    for vid, orbit in enumerate(dm.vertex_orbits()):
        for d in orbit:
            dart_to_vertex[d] = vid
    return dart_to_vertex


def _build_dart_to_face(dm: DartMap) -> list[int]:
    """Build a dart-index → face-ID lookup array.

    Parameters
    ----------
    dm : DartMap
        Input dart map.

    Returns
    -------
    list[int]
        Array of length ``dm.num_darts`` mapping each dart to its face ID.
    """
    dart_to_face: list[int] = [-1] * dm.num_darts
    for fid, orbit in enumerate(dm.face_orbits()):
        for d in orbit:
            dart_to_face[d] = fid
    return dart_to_face


class PlanarEmbeddingView:
    """Thin view over a DartMap exposing integer-indexed embedding queries.

    Parameters
    ----------
    dm : DartMap
        The combinatorial map to wrap.  Must represent a connected closed
        surface; planarity is not validated here.

    Notes
    -----
    All vertex and face identifiers are integers ``0..V-1`` / ``0..F-1``
    assigned in ``vertex_orbits()`` / ``face_orbits()`` scan order.

    Construction takes ``O(n)`` time and space where ``n = dm.num_darts``.
    """

    def __init__(self, dm: DartMap) -> None:
        """Initialise the view by building dart-to-vertex/face lookups.

        Parameters
        ----------
        dm : DartMap
            Input dart map.
        """
        self._dm = dm

        vertex_orbits = dm.vertex_orbits()
        self._vertex_to_rep: list[int] = [orbit[0] for orbit in vertex_orbits]
        self._dart_to_vertex: list[int] = _build_dart_to_vertex(dm)

        face_orbits = dm.face_orbits()
        self._face_to_rep: list[int] = [orbit[0] for orbit in face_orbits]
        self._dart_to_face: list[int] = _build_dart_to_face(dm)

    # ------------------------------------------------------------------
    # Counts
    # ------------------------------------------------------------------

    @property
    def num_vertices(self) -> int:
        """Number of vertices in the dart map.

        Returns
        -------
        int
            Count of vertex orbits.
        """
        return len(self._vertex_to_rep)

    @property
    def num_faces(self) -> int:
        """Number of faces in the dart map.

        Returns
        -------
        int
            Count of face orbits.
        """
        return len(self._face_to_rep)

    @property
    def num_edges(self) -> int:
        """Number of edges in the dart map.

        Returns
        -------
        int
            Equal to ``dm.num_edges``.
        """
        return self._dm.num_edges

    # ------------------------------------------------------------------
    # Vertex queries
    # ------------------------------------------------------------------

    def ordered_neighbors(self, v: int) -> list[int]:
        """Return neighbours of ``v`` in CCW cyclic order.

        Parameters
        ----------
        v : int
            Vertex ID in ``[0, num_vertices)``.

        Returns
        -------
        list[int]
            Neighbour vertex IDs visited in the ``sigma`` orbit order of
            ``v``'s representative dart, which is CCW in the plane.

        Raises
        ------
        IndexError
            If ``v`` is out of range.
        """
        if not 0 <= v < self.num_vertices:
            raise IndexError(
                f"Vertex ID {v} out of range [0, {self.num_vertices})."
            )
        rep = self._vertex_to_rep[v]
        return [
            self._dart_to_vertex[self._dm.alpha[d]]
            for d in vertex_darts(self._dm, rep)
        ]

    def degree(self, v: int) -> int:
        """Return the number of edges incident to vertex ``v``.

        Parameters
        ----------
        v : int
            Vertex ID in ``[0, num_vertices)``.

        Returns
        -------
        int
            Degree of ``v``.

        Raises
        ------
        IndexError
            If ``v`` is out of range.
        """
        if not 0 <= v < self.num_vertices:
            raise IndexError(
                f"Vertex ID {v} out of range [0, {self.num_vertices})."
            )
        return sum(1 for _ in vertex_darts(self._dm, self._vertex_to_rep[v]))

    # ------------------------------------------------------------------
    # Face queries
    # ------------------------------------------------------------------

    def face_boundary_vertices(self, f: int) -> list[int]:
        """Return boundary vertices of face ``f`` in phi orbit order.

        Parameters
        ----------
        f : int
            Face ID in ``[0, num_faces)``.

        Returns
        -------
        list[int]
            Vertex IDs visited in the ``phi`` orbit order of ``f``'s
            representative dart.  For the outer face (listed CW in
            ``from_face_lists``), this order is **CW** in the plane.  For
            interior faces (listed CCW), this order is **CCW**.

        Raises
        ------
        IndexError
            If ``f`` is out of range.
        """
        if not 0 <= f < self.num_faces:
            raise IndexError(
                f"Face ID {f} out of range [0, {self.num_faces})."
            )
        rep = self._face_to_rep[f]
        return [
            self._dart_to_vertex[d] for d in face_darts(self._dm, rep)
        ]

    # ------------------------------------------------------------------
    # Dart-level access (used by Phase 8 canonical ordering)
    # ------------------------------------------------------------------

    def vertex_of_dart(self, d: int) -> int:
        """Return the vertex ID of the source vertex of dart ``d``.

        Parameters
        ----------
        d : int
            Dart index in ``[0, dm.num_darts)``.

        Returns
        -------
        int
            Vertex ID of the vertex whose ``sigma`` orbit contains ``d``.

        Raises
        ------
        IndexError
            If ``d`` is out of range.
        """
        self._dm.validate_dart(d)
        return self._dart_to_vertex[d]

    def face_of_dart(self, d: int) -> int:
        """Return the face ID of the face whose ``phi`` orbit contains ``d``.

        Parameters
        ----------
        d : int
            Dart index in ``[0, dm.num_darts)``.

        Returns
        -------
        int
            Face ID.

        Raises
        ------
        IndexError
            If ``d`` is out of range.
        """
        self._dm.validate_dart(d)
        return self._dart_to_face[d]
