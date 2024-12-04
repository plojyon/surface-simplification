using DataStructures
using FileIO
using Meshes
using MeshIO
using Makie
using GLMakie

# buni = load("objects/bunny/reconstruction/bun_zipper.ply")
# Makie.mesh(buni)

function isSafe(a, b)
    return true
end

# errorfunction = (a, b) -> a < b  # Min-heap
# ordering = FunctorOrdering(errorfunction)

struct SimplicialComplex2D
    coords::Array{Vector{Float64}}
    vertices::Set{Int}
    edges::Set{Tuple{Int,Int}}
    triangles::Set{Tuple{Int,Int,Int}}

    _vertex_to_edges::Dict{Int,Vector{Int}}
    _edge_to_triangles::Dict{Int,Vector{Int}}

    function SimplicialComplex2D(triangles::Vector{Tuple{Float64,Float64,Float64}})
        coords = Vector{Vector{Float64}}()
        vertices = Vector{Vector{Int}}()
        edges = Vector{Tuple{Int,Int}}()
        triangles = Vector{Tuple{Int,Int,Int}}()

        _vertex_to_edges = Dict{Int,Vector{Int}}()
        _edge_to_triangles = Dict{Int,Vector{Int}}()

        for triangle in triangles
            v1, v2, v3 = triangle
            push!(coords, v1)
            push!(coords, v2)
            push!(coords, v3)

            push!(vertices, [length(coords) - 2, length(coords) - 1, length(coords)])
            push!(edges, (length(coords) - 2, length(coords) - 1))
            push!(edges, (length(coords) - 1, length(coords)))
            push!(edges, (length(coords) - 2, length(coords)))

            push!(triangles, (length(coords) - 2, length(coords) - 1, length(coords)))

            for edge in [(length(coords) - 2, length(coords) - 1), (length(coords) - 1, length(coords), length(coords) - 2)]
                if haskey(_edge_to_triangles, edge)
                    push!(_edge_to_triangles[edge], length(triangles))
                else
                    _edge_to_triangles[edge] = [length(triangles)]
                end
            end

            for vertex in [length(coords) - 2, length(coords) - 1, length(coords)]
                if haskey(_vertex_to_edges, vertex)
                    push!(_vertex_to_edges[vertex], length(edges))
                else
                    _vertex_to_edges[vertex] = [length(edges)]
                end
            end
        end

        new(coords, vertices, edges, triangles, _vertex_to_edges, _edge_to_triangles)
    end
end

struct Plane
    normal::Vector{Float64}
    offset::Float64
    u::Vector{Float64}

    function Plane(normal::Vector{Float64}, offset::Float64)
        new(normal, offset, [normal; offset])
    end
end

function fundamental_quadratic(planes::Vector{Plane})
    Q = Matrix{Float64}(undef, 4, 4)
    for i in 1:length(planes)
        Q += planes[i].u * planes[i].u'
    end
    return Q
end

function distance(point::Vector{Float64}, Q::Matrix{Float64})
    return point' * Q * point
end

function min_distance_point(point::Vector{Float64}, Q::Matrix{Float64})
    sol = Q \ [point; 1.0]
    return sol[1:3]
end
