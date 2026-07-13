module LatticeDatastructures

export Plaquettes, Dart
export num_plaquettes, num_neighbors
export sibling_dart, sibling_dart_ids
include("plaquettes.jl")

export EdgeType, UNMODIFIED, TAIL, CW_OF_TAIL, CCW_OF_TAIL
export EdgeLocation, TailedLattice
export num_edges, edge_location, edge_ids
include("tailedlattice.jl")

export rectangular_honeycomb_grid_on_sphere
export rectangular_honeycomb_grid_on_torus
include("recipes.jl")

end # module LatticeDatastructures
