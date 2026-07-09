using Documenter, LatticeDatastructures

makedocs(
    sitename="Lattice Datastructures",
    modules=[LatticeDatastructures],
    format=Documenter.HTML(edit_link="main"),
    pages=[
        "Manual" => "index.md",
        "Testing" => "testing.md",
    ]
)
