using Documenter, LatticeDatastructures

makedocs(
    sitename="Lattice Datastructures",
    modules=[LatticeDatastructures],
    format=Documenter.HTML(edit_link="main"),
    pages=[
        "Overview" => "index.md",
        "Recipes" => "recipes.md",
        "Visualization" => "visualization.md",
    ]
)
