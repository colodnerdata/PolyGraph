from polygraph.structures.dart_map import DartMap


def prism(n: int) -> DartMap:
    """Return the n-gonal prism with even-indexed top vertices and odd-indexed bottom vertices."""
    if n < 3:
        raise ValueError("prism(n) requires n >= 3")  # A prism must have at least a triangular base

    top = list(range(0, 2 * n, 2))        # Even indices: 0,2,4,... → vertices on the top cap
    bottom = list(range(1, 2 * n, 2))     # Odd indices: 1,3,5,... → vertices on the bottom cap

    faces = [top, bottom[::-1]]           # Two caps; reverse bottom so outward orientations match

    faces += [
        [top[i], top[(i + 1) % n], bottom[(i + 1) % n], bottom[i]]
        # Side quad connecting the ith edge of the top cap to the ith edge of the bottom cap
        # (i+1)%n wraps around to close the cycle
        for i in range(n)
    ]

    return DartMap.from_face_lists(faces, num_vertices=2 * n)  # Build the combinatorial map

from polygraph.structures.dart_map import DartMap


def antiprism(n: int) -> DartMap:
    """Return the n-gonal antiprism with even-indexed top vertices and odd-indexed bottom vertices."""
    if n < 3:
        raise ValueError("antiprism(n) requires n >= 3")  # Antiprism family starts at n = 3

    # Vertex numbering convention:
    # even indices → top cap
    # odd indices  → bottom cap
    #
    # Example for n=5:
    # top    = [0,2,4,6,8]
    # bottom = [1,3,5,7,9]
    top = list(range(0, 2 * n, 2))
    bottom = list(range(1, 2 * n, 2))

    # Add the two n-gon caps.
    #
    # The bottom face is reversed so that adjacent faces traverse
    # shared edges in opposite directions, which is typically
    # required by combinatorial map constructors.
    faces = [top, bottom[::-1]]

    # First set of side triangles.
    #
    # Each triangle connects a top edge to the bottom vertex beneath it.
    #
    # Pattern:
    #   [top_i, top_{i+1}, bottom_i]
    #
    # The modulo wraps around to close the cycle.
    faces += [
        [top[i], top[(i + 1) % n], bottom[i]]
        for i in range(n)
    ]

    # Second set of side triangles.
    #
    # These complete the "twisted band" of triangles around the polyhedron.
    #
    # Pattern:
    #   [top_{i+1}, bottom_{i+1}, bottom_i]
    faces += [
        [top[(i + 1) % n], bottom[(i + 1) % n], bottom[i]]
        for i in range(n)
    ]

    # Build the DartMap from the face list.
    # Total vertices = 2n
    return DartMap.from_face_lists(faces, num_vertices=2 * n)
