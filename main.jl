using DataStructures
using FileIO
using Meshes
using MeshIO
using Makie
using GLMakie
using GeometryBasics

include("contractions.jl")
include("visualisation.jl")


function CreateSimplicialComplex2D(mesh)
    coords = Vector{Vector{Float32}}()
    anticoords = Dict{Vector{Float32},Int}()

    triangles_ = Vector{Tuple{Int,Int,Int}}()

    for triangle in mesh
        for vertex in triangle
            if !haskey(anticoords, vertex)
                push!(coords, collect(vertex))
                anticoords[collect(vertex)] = length(coords)
            end
        end
        push!(triangles_, (anticoords[collect(triangle[1])], anticoords[collect(triangle[2])], anticoords[collect(triangle[3])]))
    end

    SimplicialComplex2D(triangles_, coords)
end

function loadsc(path)
    data = load(path)
    data_triangles = Array([Array([Tuple([z for z in y]) for y in x]) for x in data])
    CreateSimplicialComplex2D(data_triangles)
end


bunidata = loadsc("bunny/reconstruction/bun_zipper.ply")
buni = initialContractedSimplicialComplex2D(bunidata)

# for i in 1:100
#     contract!(buni)
# end

visualize(buni)
