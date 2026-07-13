"""
    EdgeType UNMODIFIED TAIL CW_OF_TAIL CCW_OF_TAIL

The four 'types' of edges in the tailed lattice. See `TailedLattice` for more.

# Values
- `UNMODIFIED`: an edge from the untailed lattice that does not host a tail; in
  other words it is not a tail and does not share a vertex with a tail
- `TAIL`: a tail
- `CW_OF_TAIL`: an edge created by the insertion of a tail, located clockwise
  around the vertex it shares with the tail that created it
- `CCW_OF_TAIL`: an edge created by the insertion of a tail, located counter-
  clockwise around the vertex it shares with the tail that created it
"""
@enum EdgeType UNMODIFIED TAIL CW_OF_TAIL CCW_OF_TAIL

"""
    EdgeLocation

Bundles the information needed to uniquely identify the location of an edge in
the tailed lattice. See `TailedLattice` for more.

# Fields
- `dart::Dart`: the side of the untailed lattice edge this edge corresponds to
- `type::EdgeType`: the relation of this edge to the side of its corresponding
  untailed lattice edge
"""
struct EdgeLocation
    dart::Dart
    type::EdgeType
end # struct EdgeLocation

"""
    TailedLattice

Plaquette- and edge-centric representation of the tailed lattice built on top of
`Plaquettes`.

We restrict each plaquette to have zero or one tails in it.

Each tail in the tailed lattice is 'attached' to one side of a host edge in the
untailed lattice. We restrict each untailed lattice edge to hosting one tail at
most, which means that each tail can be identified with a `Dart` representing
the side of the host edge it is attached to.

Furthermore, attaching a tail to a host edge splits the host into two edges, one
clockwise of and one counterclockwise of the tail, going around the vertex they
share with the tail. Thus the introduction of a tail creates three edges from
one, and every edge in the tailed lattice can be given one of four `EdgeType`s
depending on its relation to the introduced tails: either it is one of the three
edges created by adding a tail, or it is an unmodified edge from the untailed
lattice which does not host a tail.

We give every edge ``E`` in the tailed lattice a globally unique edge id (as
opposed to the plaquette-unique dart ids in the untailed lattice), which can be
used to sample from the set of lattice edges. From the above, we see that the
location of ``E`` can be uniquely specified by the combination of a `Dart` and
an `EdgeType`: the `Dart` identifies a side of an edge ``F`` of the untailed
lattice, and the `EdgeType` identifies either that ``E`` is ``F`` unmodified, or
that ``E`` is one of the derivative three edges created by inserting a tail on
the selected side of ``F``. These two pieces of information are bundled together
into an `EdgeLocation`.

In the unmodified case, the choice of which `Dart` to use to represent ``F`` is
arbitrary, but in the modified case swapping the `Dart` will swap which
plaquette we interpret the tail to be 'going into', and which edge is clockwise
vs counterclockwise of it.

# Fields
- `_plaquettes::Plaquettes`: represents the plaquette structure of the lattice
- `_tails::Dict{Int, Int}`: a collection of `Dart`s corresponding to tails in
  the lattice, represented as a mapping from `plaq_id` -> `dart_id` to allow
  `O(1)` lookups
- `_edges::Dict{Int, EdgeLocation}`: a mapping from a global edge_id identifying
  one specific edge to an `EdgeLocation` identifying its location relative to
  the untailed lattice
"""
struct TailedLattice
    _plaquettes::Plaquettes
    _tails::Dict{Int,Int}
    _edges::Dict{Int,EdgeLocation}

    function TailedLattice(p::Plaquettes, tails::Dict{Int,Int})
        # validate tails
        for (plaq_id, dart_id) in tails
            1 <= plaq_id <= num_plaquettes(p) ||
                throw(ArgumentError("tail in nonexistent plaquette $plaq_id"))
            1 <= dart_id <= num_neighbors(p, plaq_id) ||
                throw(ArgumentError("tail on nonexistent dart $dart_id of plaquette $plaq_id"))
            sib = sibling_dart(p, plaq_id, dart_id)
            haskey(tails, sib.plaq_id) && tails[sib.plaq_id] == sib.dart_id &&
                throw(ArgumentError("sibling darts ($plaq_id, $dart_id) / ($(sib.plaq_id), $(sib.dart_id)) both have tails"))
        end
        # assign edges
        next_edge_id = 1
        edges = Dict{Int,EdgeLocation}()
        visited = Set{Dart}()
        for plaq_id in 1:num_plaquettes(p)
            for dart_id in 1:num_neighbors(p, plaq_id)
                # visit each edge once, by processing both darts at once
                dart = Dart(plaq_id, dart_id)
                dart ∈ visited && continue
                sib_dart = sibling_dart(p, dart.plaq_id, dart.dart_id)
                push!(visited, dart)
                push!(visited, sib_dart)
                # determine if there's a tail on either side of this edge
                if get(tails, plaq_id, nothing) == dart_id
                    taildart = dart
                elseif get(tails, sib_dart.plaq_id, nothing) == sib_dart.dart_id
                    taildart = sib_dart
                else
                    taildart = nothing
                end
                # assign edge id(s)
                if taildart === nothing
                    # register one unmodified edge
                    edges[next_edge_id] = EdgeLocation(dart, UNMODIFIED)
                    next_edge_id += 1
                else
                    # register three modified edges
                    edges[next_edge_id] = EdgeLocation(taildart, CCW_OF_TAIL)
                    next_edge_id += 1
                    edges[next_edge_id] = EdgeLocation(taildart, TAIL)
                    next_edge_id += 1
                    edges[next_edge_id] = EdgeLocation(taildart, CW_OF_TAIL)
                    next_edge_id += 1
                end
            end
        end

        new(p, tails, edges)
    end # constructor TailedLattice
end # struct TailedLattice

"""
    TailedLattice(p::Plaquettes, tails::Dict{Int, Int})

Inner constructor which accepts and validates `tails` directly.

Because `tails` is a dictionary, each plaquette can only have zero or one tails.

Two additional properties are verified:
- all tails are on valid edges in valid plaquettes
- sibling darts don't both have tails
"""
TailedLattice(p::Plaquettes, tails::Dict{Int, Int})

"""
    num_edges(tl::TailedLattice)

Return the number of edges in `tl`.
"""
num_edges(tl::TailedLattice) = length(tl._edges)

"""
    edge_location(tl::TailedLattice, edge_id::Int)

Return the `EdgeLocation` for edge `edge_id` in `tl.`
"""
edge_location(tl::TailedLattice, edge_id::Int) = tl._edges[edge_id]

"""
    edge_ids(tl::TailedLattice)

Return an iterator over the ids of the edges in `tl`.
"""
edge_ids(tl::TailedLattice) = 1:num_edges(tl)
