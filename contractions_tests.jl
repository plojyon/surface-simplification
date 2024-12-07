include("contractions.jl")

function print_struct_properties(obj::T) where T
    println("Properties of $(typeof(obj)):")
    for field in fieldnames(T)
        value = getfield(obj, field)
        println("  $field: $value")
    end
end

coords = [Vector([rand()*100, rand()*100,rand()*100]) for i in 1:9]
torus = [(2,1,4), (4,1,3), (5,4,7), (6,5,8), (7,4,6), (7,1,8), (8,2,9), (7,9,1), (1,2,8), (4,5,2), (7,8,5), (3,6,4), (6,9,7), (1,9,3), (8,9,6), (3,2,5), (5,6,3), (2,3,9)]

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

println(simplicial_complex)
contracted = initialContractedSimplicialComplex2D(simplicial_complex)
contracted = contract(contracted, (1, 2))
println(contracted.contracted)