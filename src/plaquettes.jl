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
"""
struct Plaquettes
    _siblings::Vector{Vector{Dart}}

    function Plaquettes(siblings::Vector{Vector{Dart}})
        for plaq_id in eachindex(siblings)
            for dart_id in eachindex(siblings[plaq_id])
                # looking up every Dart in siblings enforces id validity
                D1 = Dart(plaq_id, dart_id)
                D2 = siblings[plaq_id][dart_id]
                # no self-neighbors check
                D1.plaq_id == D2.plaq_id &&
                    throw(ArgumentError("sibling darts $D1 and $D2 are in the same plaquette"))
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
plaquette must not neighbor itself.

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

Return tuples of `dart_ids`s which correspond to edges shared between
`plaq1` and `plaq2` in `p`.
"""
function sibling_dart_ids(p::Plaquettes, plaq_id1::Int, plaq_id2::Int)
    pairs = Vector{Tuple{Int, Int}}()
    for (dart_id1, dart2) in enumerate(p._siblings[plaq_id1])
        dart2.plaq_id == plaq_id2 && push!(pairs, (dart_id1, dart2.dart_id))
    end
    pairs
end # function sibling_dart_ids
