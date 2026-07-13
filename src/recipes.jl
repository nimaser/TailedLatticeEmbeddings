"""
    rectangular_honeycomb_grid_on_sphere()

Create a `TailedLattice` whose topology is a `rows`x`cols` grid of hexagons,
indexed as `(r, c)` starting from the top left, where all hexagons are oriented
so that two of their sides are 'vertical' and the tails are on the left vertical
side of each plaquette. The rows would have alternating starting points at 0 and
0.5 hexagon widths from the left, if laid out in a regular grid, or in other
words the grid topology is as rectangular as possible rather than being a non-
rectangular parallelogram, where every row would start 0.5 hexagon widths to the
right of the row above.

The plaquettes on the upper, lower, left and right boundaries of the grid all
border a single plaquette which occupies the 'backside' of the sphere and has
`(2*cols) + (2*cols) + (2*(rows) - 1) + (2*(rows) - 1) = 4*cols + 4*rows - 2`
sides. Thus the lattice will have a total of `rows*cols + 1` plaquettes.
"""
function rectangular_honeycomb_grid_on_sphere(rows::Int, cols::Int)

end

"""
    rectangular_honeycomb_grid_on_torus()

Create a `TailedLattice` whose topology is a `rows`x`cols` grid of hexagons,
indexed as `(r, c)` starting from the top left, where all hexagons are oriented
so that two of their sides are 'vertical' and the tails are on the left vertical
side of each plaquette. The rows would have alternating starting points at 0 and
0.5 hexagon widths from the left, if laid out in a regular grid, or in other
words the grid topology is as rectangular as possible rather than being a non-
rectangular parallelogram, where every row would start 0.5 hexagon widths to the
right of the row above.

In order to make the toroidal boundary condition work nicely, we require that
the lattice have an even number of rows. In this way, the top row alternates
with the bottom row, meaning the pointy bits of the top one slot between the
pointy bits of the bottom one. Due to the rectangular grid, the left and right
columns also line up. Thus the lattice can be imagined as taking the grid and
rolling it up along the top/bottom and left/right boundaries, as one would
do to make a torus out of a rectangular sheet of paper. This lattice has a total
of `rows*cols` plaquettes.
"""
function rectangular_honeycomb_grid_on_torus(rows::Int, cols::Int)

end
