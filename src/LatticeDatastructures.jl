module LatticeDatastructures

export Plaquettes, Dart
export num_plaquettes, num_neighbors
export sibling_dart, sibling_dart_ids

export EdgeType, UNMODIFIED, TAIL, CW_OF_TAIL, CCW_OF_TAIL
export EdgeLocation, TailedLattice
export num_edges, edge_location, edge_ids

################################################################################
# Plaquettes
################################################################################

"""
    Dart

Uniquely identifies one 'side' of an edge shared between two plaquettes in a
`Plaquettes`, specifically the `plaq_id` side of the edge with `plaq_id`- local
label `edge_id`. See `Plaquettes` for more.

# Fields
- `plaq_id::Int`: must be positive and unique within each `Plaquettes`
- `dart_id::Int`: must be positive and unique within each plaquette
"""
struct Dart
    plaq_id::Int
    dart_id::Int
    Dart(plaq_id::Int, dart_id::Int) =
        plaq_id > 0 && dart_id > 0 ?
        new(plaq_id, dart_id) :
        throw(ArgumentError("plaq and dart ids must be positive"))
end # struct Dart

"""
    Plaquettes

Combinatorial map of the dual of the embedding of the untailed lattice. This
representation facilitates querying the topological arrangement of the neighbors
of each plaquette.

We assign each plaquette in the lattice a global id, and assign each of its
edges a plaquette-local id. With this scheme, each edge has two plaquette-local
ids, one for each plaquette it is shared by. The combination of a plaquette id
and a plaquette-local edge id defines one 'side' of an edge, and is bundled
together into a `Dart`. To avoid confusion between globally- and locally-unique
edge ids, we call the plaquette-local edge ids 'dart ids'.

Therefore, an edge being shared between two plaquettes can be indicated by
pairing together the two `Dart`s associated with the two sides of the shared
edge. These two are termed "sibling" `Dart`s.

By storing a plaquette's `Dart`s in clockwise order, the local topology of the
lattice around that plaquette can be stored. Doing this for all plaquettes
stores the topology of the entire lattice.

We can simultaneously and simply achieve both of these objectives with one
datastructure:
- first, we assign plaquette ids from 1 to the number of lattice plaquettes
- next, we assign dart ids for plaquette P in clockwise order from 1 to the
  number of edges of P
- This allows us to use `plaq_id` and `dart_id` as indices: we store the sibling
  of `Dart(plaq_id, dart_id)` at `_siblings[plaq_id][dart_id]`
- Furthermore, because we can get the previous/next `Dart` by subtracting/adding
    1, mod the number of edges of P, to the current `Dart`'s `dart_id`, storing the
    darts in an array gives us the ordering of neighboring plaquettes for free

# Fields
- `_siblings::Vector{Vector{Dart}}`: a mapping between clockwise-ordered `Darts`
  in a plaquette and their siblings in neighboring plaquettes

...
"""
struct Plaquettes
    _siblings::Vector{Vector{Dart}}

    function Plaquettes(siblings::Vector{Vector{Dart}})
        for plaq_id in eachindex(siblings)
            neighbors = Set{Int}()
            for dart_id in eachindex(siblings[plaq_id])
                # looking up every Dart in siblings enforces id validity
                D1 = Dart(plaq_id, dart_id)
                D2 = siblings[plaq_id][dart_id]
                # no self-neighbors check
                D1.plaq_id == D2.plaq_id &&
                    throw(ArgumentError("sibling darts $D1 and $D2 are in the same plaquette"))
                # shared edge uniqueness check
                D2.plaq_id ∈ neighbors &&
                    throw(ArgumentError("$(D2.plaq_id) neighbors $plaq_id multiple times"))
                push!(neighbors, D2.plaq_id)
                # bijectivity check
                D2_sib = siblings[D2.plaq_id][D2.dart_id]
                D1 == D2_sib ||
                    throw(ArgumentError("sibling mapping is not bijective between $D1 and $D2"))
            end
        end
        new(siblings)
    end # constructor Plaquettes

end # struct Plaquettes

# inner constructor docstring
"""
    Plaquettes(siblings::Vector{Vector{Dart}})

Inner constructor which accepts and validates `siblings` directly.

Several properties are verified:
- id validity: plaquettes have ids from 1 to the number of plaquettes, and
  darts have labels from 1 to the number of darts in the their plaquete
- no self-neighbors: sibling darts cannot belong to the same plaquette
- shared edge uniqueness: no pair of plaquettes can share more than 1 edges
- dart uniqueness: each dart should appear in `siblings` exactly one time,
  when it is mapped to by its sibling
- bijectivity: if dart D1 maps to dart D2, then dart D2 also maps to dart D1
    - note that dart uniqueness is implied by bijectivity
"""
Plaquettes(siblings::Vector{Vector{Dart}})

"""
    Plaquettes(neighbors::Vector{Vector{Int}})

An alternative constructor which accepts a specification of the neighbors of
each plaquette, converts it to a siblings array, and passes it to the inner
constructor.

Performs a check that in each pair of plaquettes, they both list the other
as a neighbor the same number of times. Leaves all other validation to the
main constructor: in particular, its implementation does not assume that a
plaquette must not neighbor itself or must only share one edge with each
neighbor.

Each neighbor of a plaquette X is assigned one dart in X: the nth listed
neighbor will be assigned `dart_id` n.

# Fields
- `neighbors::Vector{Vector{Int}}`: `neighbors[plaq_id]` should be a
  clockwise-ordered list of the plaquette ids of the neighbors of `plaq_id`
"""
function Plaquettes(neighbors::Vector{Vector{Int}})
    # initialize siblings
    num_plaqs = length(neighbors)
    siblings = [Vector{Dart}(undef, length(neighbor_list)) for neighbor_list in neighbors]
    # iterate over unordered pairs of plaquettes {x, y}
    for plaq_x in 1:num_plaqs, plaq_y in (plaq_x):num_plaqs
        # find where x and y list each other as neighbors
        inds_y_in_x = findall(==(plaq_y), neighbors[plaq_x])
        inds_x_in_y = findall(==(plaq_x), neighbors[plaq_y])
        length(inds_y_in_x) == length(inds_x_in_y) ||
            throw(ArgumentError("plaquettes $plaq_x and $plaq_y list each other as neighbors different numbers of times"))
        # assign matching pairs of darts: reverse the latter so traversal
        # direction along the shared sides is the same: draw to understand
        for (dart_id_x, dart_id_y) in zip(inds_y_in_x, reverse(inds_x_in_y))
            siblings[plaq_x][dart_id_x] = Dart(plaq_y, dart_id_y)
            siblings[plaq_y][dart_id_y] = Dart(plaq_x, dart_id_x)
        end
    end
    Plaquettes(siblings)
end # constructor Plaquettes

"""
    num_plaquettes(p::Plaquettes)

Return the number of plaquettes in `p`.
"""
num_plaquettes(p::Plaquettes) = length(p._siblings)

"""
    num_neighbors(p::Plaquettes, plaq_id::Int)

Return the number of neighboring plaquettes to `plaq_id` in `p`.
"""
num_neighbors(p::Plaquettes, plaq_id::Int) = length(p._siblings[plaq_id])

"""
    sibling_dart(p::Plaquettes, plaq_id::Int, dart::Int)

Return the sibling to the `dart_id`th dart in `plaq_id` in `p`.
"""
sibling_dart(p::Plaquettes, plaq_id::Int, dart_id::Int) = p._siblings[plaq_id][dart_id]

"""
    sibling_dart_ids(p::Plaquettes, plaq_id1::Int, plaq_id2::Int)

Return the `dart_id`s of the pair of darts corresponding to the
edge shared between `plaq1` and `plaq2` in `p`.
"""
function sibling_dart_ids(p::Plaquettes, plaq_id1::Int, plaq_id2::Int)
    for (dart_id1, dart2) in enumerate(p._siblings[plaq_id1])
        dart2.plaq_id == plaq_id2 && return dart_id1, dart2.dart_id
    end
    nothing
end # function sibling_dart_ids

################################################################################
# TailedLattice
################################################################################

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
  the lattice, represented as a mapping from plaq_id -> dart_id to allow O(1)
  lookups
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
                    continue
                end
                # assign edge id(s)
                if t === nothing
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
num_edges(tl::TailedLattice) = length(t._edges)

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

end # module LatticeDatastructures
