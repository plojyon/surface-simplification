include("contractions.jl")

function print_struct_properties(obj::T) where {T}
    println("Properties of $(typeof(obj)):")
    for field in fieldnames(T)
        value = getfield(obj, field)
        println("  $field: $value")
    end
end

coords::Vector{Vector{Float32}} = [Vector([rand() * 100, rand() * 100, rand() * 100]) for i in 1:9]
torus = [(2, 1, 4), (4, 1, 3), (5, 4, 7), (6, 5, 8), (7, 4, 6), (7, 1, 8), (8, 2, 9), (7, 9, 1), (1, 2, 8), (4, 5, 2), (7, 8, 5), (3, 6, 4), (6, 9, 7), (1, 9, 3), (8, 9, 6), (3, 2, 5), (5, 6, 3), (2, 3, 9)]

simplicial_complex = SimplicialComplex2D(torus, coords)
# println("Edges that contain 1")
# arr = edges(simplicial_complex)
# for edge in sort(collect(arr), by = x -> x[2])
#     if edge[1] == 1 || edge[2] == 1
#         println(edge)
#     end
# end

# println("Edges that do not contain 1")
# for edge in sort(collect(arr), by = x -> x[2])
#     if edge[1] != 1 && edge[2] != 1
#         println(edge)
#     end
# end

# println("Triangles that contain 1")
# arr = triangles(simplicial_complex)
# for triangle in sort(collect(arr), by = x -> x[2])
#     if triangle[1] == 1 || triangle[2] == 1 || triangle[3] == 1
#         println(triangle)
#     end
# end

# println("Triangles that do not contain 1")
# for triangle in sort(collect(arr), by = x -> x[2])
#     if triangle[1] != 1 && triangle[2] != 1 && triangle[3] != 1
#         println(triangle)
#     end
# end

# println("\nReplacing vertex 1 with vertex 10\n")
# contracted = initialContractedSimplicialComplex2D(simplicial_complex)
# replaceVertexInEdgesAndTriangles(contracted.contracted, 1, 10)

# println("Edges that contain 10:")
# arr = edges(contracted.contracted)
# for edge in sort(collect(arr), by = x -> x[1])
#     if edge[1] == 10 || edge[2] == 10
#         println(edge)
#     end
# end

# println("Edges that do not contain 10:")
# for edge in sort(collect(arr), by = x -> x[2])
#     if edge[1] != 10 && edge[2] != 10
#         println(edge)
#     end
# end

# println("Triangles that contain 10:")
# arr = triangles(contracted.contracted)
# for triangle in sort(collect(arr), by = x -> x[1])
#     if triangle[1] == 10 || triangle[2] == 10 || triangle[3] == 10
#         println(triangle)
#     end
# end

# println("Triangles that do not contain 10:")
# for triangle in sort(collect(arr), by = x -> x[2])
#     if triangle[1] != 10 && triangle[2] != 10 && triangle[3] != 10
#         println(triangle)
#     end
# end

# println(" === contract tests === ")
# println(simplicial_complex)
# contracted = initialContractedSimplicialComplex2D(simplicial_complex)
# c = contract!(contracted, (1, 2))
# println(" - the edge (1, 2) does not exist anymore")
# println("result: ", !haskey(contracted.contracted._edge_to_triangles, (1, 2)))

# println(" - all edges have Q matrix")
# is_missing = false
# for edge in edges(contracted.contracted)
#     if !haskey(contracted._contracted_edge_Q, edge)
#         global is_missing
#         println("edge $edge does not have Q matrix")
#         is_missing = true
#     end
# end
# println("result: ", !is_missing)

# println(" - all triangles have Q matrix")
# is_missing = false
# for triangle in triangles(contracted.contracted)
#     if !haskey(contracted._contracted_triangle_Q, triangle)
#         global is_missing
#         println("triangle $triangle does not have Q matrix")
#         is_missing = true
#     end
# end
# println("result: ", !is_missing)

# println(" === contract tests === ")
# println(simplicial_complex)
# contracted = initialContractedSimplicialComplex2D(simplicial_complex)
# contract!(contracted, (1, 2))
# println(contracted.contracted)

# println(" === getlink tests === ")
# println(getlink(simplicial_complex, 1))
# println(getlink(simplicial_complex, 2))
# println(getlink(simplicial_complex, (1, 2)))

#println(error(contracted, (1, 2)))
a = 5
b = 6
x = 2
y = 9
poper = [
    (1, x, a),
    (1, x, 3),
    (x, a, b),
    (x, b, 3),
    (1, a, 4),
    (a, 4, 8),
    (a, 8, y),
    (b, 3, 7),
    (b, 7, 10),
    (b, y, 8),
    (b, y, 10),
    (7, 10, y),
    (4, 8, y)
]
poper_coordinates = [
    [Float32(0.0), Float32(1.0), Float32(0.0)],  # "1"
    [Float32(1.0), Float32(1.5), Float32(0.0)],  # "x"
    [Float32(2.0), Float32(1.0), Float32(0.0)],  # "3"
    [Float32(0.0), Float32(0.0), Float32(0.0)],  # "4"
    [Float32(1.0), Float32(0.5), Float32(0.0)],  # "a"
    [Float32(2.0), Float32(0.5), Float32(0.0)],  # "b"
    [Float32(2.0), Float32(-0.5), Float32(0.0)], # "7"
    [Float32(1.0), Float32(-1.0), Float32(0.0)], # "8"
    [Float32(0.0), Float32(-0.5), Float32(0.0)], # "y"
    [Float32(2.5), Float32(-1.5), Float32(0.0)]  # "0"
]

popercomplex = SimplicialComplex2D(poper, poper_coordinates)
println(" === poper tests === ")
contracted = initialContractedSimplicialComplex2D(popercomplex)
contract!(contracted, (a, b))
println(contracted.contracted)
