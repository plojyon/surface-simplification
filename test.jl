include("contractions.jl")
include("main.jl")

println(" - edge has 2 adjacent triangles")
is_missing = false
for edge in edges(buni.original)
    if length(buni.original._edge_to_triangles[edge]) != 2
        global is_missing
        println("edge $edge does not have 2 adjacent triangles")
        is_missing = true
    end
end
println("result: ", !is_missing)