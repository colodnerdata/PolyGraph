# Contributing

PolyGraph is a research-oriented Python library for combinatorial topology, planar graph drawing, and polyhedral realization. Contributions should favor clarity, mathematical correctness, and maintainability over cleverness.

General principles

Write the simplest correct version first. Prefer explicit code over abstraction when the abstraction makes the mathematics harder to follow.

Keep the architectural layers separate:

core algebra and permutations

combinatorial structures

traversal and topology

algorithms

geometry

visualization and export


Do not mix topology and geometry unless a module is explicitly responsible for bridging them.

Code style

Use Python 3.11+ features when they improve readability.

Follow PEP 8, with a maximum line length of 79 characters.

Use descriptive names. In mathematical code, short names are acceptable when they are standard and local, such as i, j, n, sigma, alpha, or phi. Outside tight mathematical contexts, prefer explicit names.

Prefer pure functions where practical. Hidden mutation should be avoided unless performance clearly requires it.

Keep functions small and single-purpose. If a function has multiple conceptual stages, split it.

Avoid premature optimization. Optimize only after correctness and tests are in place.

Types and signatures

Use type hints on all public functions and methods.

Annotate return types explicitly.

Prefer concrete, readable types in public APIs. Use protocols or abstract types only when they materially improve reuse.

Raise clear exceptions with layer-appropriate language. Low-level modules should use generic terms like “index” or “permutation,” not higher-level terms like “dart,” unless the module is specifically about dart maps.

Documentation style

PolyGraph uses NumPy-style docstrings.

Public functions, methods, classes, and modules should have docstrings unless they are obviously internal and trivial.

Required docstring sections

Use these when applicable:

summary line

parameters

returns or yields

raises

notes


Optional docstring sections

These are optional and should be included only when they add value:

examples

references

see also


Do not document trivial examples just to fill space. Examples should be included only when they clarify non-obvious behavior, mathematical conventions, edge cases, or expected usage.

References are encouraged for algorithms, formulas, and nonstandard terminology, but they are not required for routine helper functions.

Docstring rules

Use imperative mood in the first sentence.

Good: Return the inverse permutation.

Bad: Returns the inverse permutation.

Do not repeat type information in prose.

Good: start : int     Starting index of the orbit.

Bad: start : int     The starting index as an integer.

Put mathematical explanation in Notes, not in the summary line.

Example docstring template

def orbit(self, start: int) -> Iterator[int]:
    """Iterate over the orbit of an index under the permutation.

    Parameters
    ----------
    start : int
        Index whose orbit will be traversed.

    Yields
    ------
    int
        Successive elements of the orbit.

    Raises
    ------
    IndexError
        If ``start`` is outside the permutation domain.

    Notes
    -----
    The orbit is defined by repeated application of the permutation until
    the starting index repeats.
    """

Module docstrings

Each substantial module should begin with a short conceptual overview.

A good module docstring explains:

what the module is for

the main mathematical objects it defines

any important conventions


For foundational modules, mention key identities or conventions when helpful.

Mathematical conventions

State important conventions explicitly and keep them consistent.

Examples:

composition order for permutations

orientation conventions for faces

indexing conventions for darts

whether generators return representatives or full orbits


If a convention is easy to misunderstand, document it in the module docstring and in the relevant method docstrings.

Testing expectations

Every nontrivial contribution should include tests.

Tests should focus first on correctness and invariants, not incidental implementation details.

Prefer small, explicit fixtures. Polyhedral generators like tetrahedron, cube, prism, and antiprism should be used heavily as canonical test cases.

For mathematical code, test invariants directly. Examples include:

permutation validity

involution properties

orbit partitioning

Euler characteristic

genus

expected vertex, edge, and face counts

symmetry group orders

planarity and convexity conditions where relevant


Research-oriented guidance

When implementing an algorithm from the literature, prefer a straightforward version first, even if asymptotically suboptimal, provided it is correct and readable.

Document the source in Notes or References when the implementation follows a paper or textbook closely.

If the implementation intentionally deviates from a reference, say so.

Organization and layering

Place code in the lowest layer that can own it cleanly.

As a rule:

generic algebra belongs in top-level foundational modules

stored combinatorial objects belong in structures

graph and topology procedures belong in algorithms

coordinate-producing code belongs in geometry

rendering helpers belong in visualization

file format writers belong in export


Avoid circular dependencies. Lower layers must not depend on higher ones.

Pull request expectations

Keep pull requests focused.

A good contribution usually includes:

the implementation

tests

docstrings

a brief note explaining the mathematical idea if it is non-obvious


For larger algorithmic additions, include a short explanation of the approach and any references used.

Final rule

Favor code that a mathematically literate reader can understand six months later.
