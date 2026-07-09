## Introduction

We need a computational representation of the tailed lattice ``T``  we use to
realize the QEC code we are simulating. In particular, we need a datastructure
capable of representing ``T``'s embedding in a closed and oriented 2D manifold.

Because ``T`` is a ribbon graph, and ribbon graphs are equivalent up to smooth
deformations of their embeddings, the embedding does not have any geometric
features we need to keep track of, only topological ones: in particular, we need
to store the adjacency relations between the lattice's plaquettes, which are
needed by the decoder and curve diagram implementations. In addition, we also
need to be able to sample edges during the noise phase of the simulation, and
know what plaquettes will be affected by noise on a particular edge.

[Combinatorial maps](https://en.wikipedia.org/wiki/Combinatorial_map) are one
candidate representation. They would completely capture all of the topological
information about the embedding, and are vertex- and edge- centric, which would
facilitate easy extension to models with e.g. multiple tails per edge/plaquette
or lattice topologies that change during the simulation. However, although
plaquette adjacency queries are technically possible, they require complicated
traversals, implementation of which will be time-consuming and error-prone.
Therefore, a plaquette- and edge-centric representation would be a better choice.

[Dual graphs](https://en.wikipedia.org/wiki/Dual_graph) effectively represent
information about the plaquettes and edges of an embedded graph, and a
combinatorial map could be used to store it without needing reference to the
primal graph. Therefore, a combinatorial map of the dual graph ``D`` of ``T``
is another good candidate. However, the tails present a slight complication:
each tail in ``T`` will cause a self-edge in ``D``. Furthermore, each tail
splits its host edge into two. In principle this is okay, but being able to
assume that
- each plaquette shares at most one edge with any other plaquette
- every vertex has three different plaquettes meeting at it (equivalently:
  no plaquette can be its own neighbor)

makes the implementation of curve diagrams much easier. Therefore, we will
instead enforce these conditions on the untailed version of ``T``, ``U``, and
use a combinatorial map ``C`` of ``U``'s dual graph to store the plaquette
adjacency.

We will store the tails as additional metadata on top of ``C``. To make
bookkeeping easier, we will assume that each plaquette has only one tail, and
that each edge can only host a tail for one of the two plaquettes it is shared
by. For the purposes of our simulation this will not make a difference.


This separated arrangement means that any code which just needs to deal with
plaquettes only needs ``C``. Only edge sampling code needs the tailed lattice.

## Untailed Lattice

```@docs
Plaquettes
Plaquettes(siblings::Vector{Vector{Dart}})
Plaquettes(neighbors::Vector{Vector{Int}})
Dart
num_plaquettes
num_neighbors
sibling_dart
sibling_dart_ids
```

## Tailed Lattice

```@docs
TailedLattice
EdgeType
EdgeInfo

```
