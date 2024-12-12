println("Julia running. Loading packages")

println("DataStructures")
using DataStructures
println("FileIO")
using FileIO
println("Meshes")
using Meshes
println("MeshIO")
using MeshIO
println("GeometryBasics")
using GeometryBasics
println("ProgressBars")
using ProgressBars
println("Serialization")
using Serialization
println("JSON3")
using JSON3
println("Serialization")
using Serialization
println("JLD2")
using JLD2
println("Makie")
using Makie
println("GLMakie")
using GLMakie

println("Including files")
include("contractions.jl")
include("visualisation.jl")
include("preprocessing.jl")

function cachesc(sc, path)
    JSON3.write("$path.json", sc)
    serialize("$path.dat", sc)
    save_object("$path.jld2", sc)
end

println("Loading bunny")
#bunidata = loadsc("bunny/reconstruction/bun_zipper.ply")
bunidata = deserialize("bunidata.dat")
buni = initialContractedSimplicialComplex2D(bunidata)
fillholes!(buni.contracted)
calculateFundamentalQuadratics!(buni)

# # mark hole borders
# marked_vertices = Set{Int}()
# for edge in getborder(buni.original)
#     push!(marked_vertices, edge[1])
#     push!(marked_vertices, edge[2])
# end
# marks::Vector{GeometryBasics.Point{3,Float32}} = [GeometryBasics.Point{3,Float32}(buni.original.coords[x]) for x in marked_vertices]
# visualize(buni, marks)

# pq = PriorityQueue()
# for (edge, triangle) in ProgressBar(buni.contracted._edge_to_triangles)
#     pq[edge] = error(buni, edge)
# end

# println("Starting contractions")
# for i in ProgressBar(1:3)
#     first = dequeue!(pq)
#     while !issafe(buni.contracted, first)
#         first = dequeue!(pq)
#     end

#     new_vertex = contract!(buni, first)
#     for edge in buni.contracted._vertex_to_edges[new_vertex]
#         pq[edge] = error(buni, edge)
#         for triangle in buni.contracted._edge_to_triangles[edge]
#             pq[triangle] = error(buni, triangle)
#         end
#     end
# end

println(" === contract tests === ")
contracted = buni

println(" - all edges have Q matrix")
is_missing = false
for edge in edges(contracted.contracted)
    if !haskey(contracted._contracted_edge_Q, edge)
        global is_missing
        println("edge $edge does not have Q matrix")
        is_missing = true
    end
end
println("result: ", !is_missing)

println(" - all triangles have Q matrix")
is_missing = false
for triangle in triangles(contracted.contracted)
    if !haskey(contracted._contracted_triangle_Q, triangle)
        global is_missing
        println("triangle $triangle does not have Q matrix")
        is_missing = true
    end
end
println("result: ", !is_missing)
