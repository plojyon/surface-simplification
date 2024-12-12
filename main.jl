println("Julia running. Loading packages")

println("DataStructures")
using DataStructures
println("FileIO")
using FileIO
println("Meshes")
using Meshes
println("MeshIO")
using MeshIO
println("Makie")
using Makie
println("GLMakie")
using GLMakie
println("GeometryBasics")
using GeometryBasics
println("ProgressBars")
using ProgressBars

println("Including files")
include("contractions.jl")
include("visualisation.jl")

function findPointInRadius(vertex, vertexSet, radius)
    for v in vertexSet
        if norm(vertex - v) < radius
            return v
        end
    end
    return nothing
end

function CreateSimplicialComplex2D(mesh)
    coords = Vector{Vector{Float32}}()
    anticoords = Dict{Vector{Float32},Int}()

    triangles_ = Vector{Tuple{Int,Int,Int}}()

    for triangle in ProgressBar(mesh)
        for vertex in triangle
            arr_vertex = collect(vertex)
            close_point = findPointInRadius(arr_vertex, coords, 0.0001)

            if close_point == nothing
                push!(coords, arr_vertex)
                anticoords[arr_vertex] = length(coords)
            else
                anticoords[arr_vertex] = anticoords[close_point]
            end
        end
        push!(triangles_, (anticoords[collect(triangle[1])], anticoords[collect(triangle[2])], anticoords[collect(triangle[3])]))
    end

    SimplicialComplex2D(triangles_, coords)
end

function loadsc(path)
    data = load(path)
    data_triangles = Array([Array([Tuple([round(z, digits=5) for z in y]) for y in x]) for x in data])
    CreateSimplicialComplex2D(data_triangles)
end


println("Loading bunny")
bunidata = loadsc("bunny/reconstruction/bun_zipper.ply")
buni = initialContractedSimplicialComplex2D(bunidata)

visualize(buni)

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

