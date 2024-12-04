using LinearAlgebra

const Triangle = Tuple{Int,Int,Int}
const Edge = Tuple{Int,Int}

function sortEdge(edge::Edge)
    return edge[1] < edge[2] ? edge : (edge[2], edge[1])
end

function sortTriangle(triangle::Triangle)
    return Tuple(sort(collect(triangle)))
end

struct SimplicialComplex2D
    coords::Array{Vector{Float64}}
    vertices::Set{Int}
    edges::Set{Edge}
    triangles::Set{Triangle}

    _vertex_to_edges::Dict{Int,Array{Edge}}
    _edge_to_triangles::Dict{Edge,Array{Triangle}}

    function SimplicialComplex2D(triangles::Array{Triangle}, coords::Array{Vector{Float64}})
        vertices = Set{Int}()
        edges = Set{Edge}()
        triangles_sorted = Set{Triangle}()

        _vertex_to_edges = Dict{Int,Array{Edge}}()
        _edge_to_triangles = Dict{Edge,Array{Triangle}}()

        for triangle in triangles
            triangle = sortTriangle(triangle)
            push!(triangles_sorted, triangle)
            for edge in [(triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1])]
                edge = sortEdge(edge)
                push!(edges, edge)
                if haskey(_edge_to_triangles, edge)
                    push!(_edge_to_triangles[edge], triangle)
                else
                    _edge_to_triangles[edge] = [triangle]
                end

                for vertex in edge
                    push!(vertices, vertex)
                    if haskey(_vertex_to_edges, vertex)
                        push!(_vertex_to_edges[vertex], edge)
                    else
                        _vertex_to_edges[vertex] = [edge]
                    end
                end
            end
        end

        new(coords, vertices, edges, triangles_sorted, _vertex_to_edges, _edge_to_triangles)
    end
end

struct ContractedSimplicialComplex2D
    original::SimplicialComplex2D
    contracted::SimplicialComplex2D
    vertex_contracted_to_original::Dict{Int,Set{Int}}

    _contracted_triangle_Q::Dict{Triangle,Matrix{Float64}}
    _contracted_edge_Q::Dict{Edge,Matrix{Float64}}
    _contracted_vertex_Q::Dict{Int,Matrix{Float64}}

    function ContractedSimplicialComplex2D(original, contracted, vertex_contracted_to_original)
        _contracted_triangle_Q = Dict{Triangle,Matrix{Float64}}()
        _contracted_edge_Q = Dict{Edge,Matrix{Float64}}()
        _contracted_vertex_Q = Dict{Int,Matrix{Float64}}()

        new(original, contracted, vertex_contracted_to_original, _contracted_triangle_Q, _contracted_edge_Q, _contracted_vertex_Q)
    end
end

function initialContractedSimplicialComplex2D(K::SimplicialComplex2D)
    vertex_contracted_to_original = Dict{Int,Set{Int}}()
    for vertex in K.vertices
        vertex_contracted_to_original[vertex] = Set([vertex])
    end
    contracted = ContractedSimplicialComplex2D(K, deepcopy(K), vertex_contracted_to_original)

    calculateFundamentalQuadratics(contracted)
    return contracted
end

# fills in the _triangle_Q, _edge_Q, and _vertex_Q fields of K
function calculateFundamentalQuadratics(K::ContractedSimplicialComplex2D)
    for triangle in K.contracted.triangles
        K._contracted_triangle_Q[triangle] = fundamental_quadratic([planeFromTriangle(triangle, K.contracted.coords)])
    end

    for edge in K.contracted.edges
        neigh_triangles = K.contracted._edge_to_triangles[edge]
        K._contracted_edge_Q[edge] = sum([K._contracted_triangle_Q[triangle] for triangle in neigh_triangles])
    end

    for vertex in K.contracted.vertices
        neigh_edges = K.contracted._vertex_to_edges[vertex]
        neigh_triangles = Set()
        for edge in neigh_edges
            neigh_triangles = neigh_triangles ∪ K.contracted._edge_to_triangles[edge]
        end
        K._contracted_vertex_Q[vertex] = sum([K._contracted_triangle_Q[triangle] for triangle in neigh_triangles])
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

function planeFromTriangle(triangle::Triangle, coords::Array{Vector{Float64}})
    #TODO [CHATGPT] confirm code

    a = coords[triangle[1]]
    b = coords[triangle[2]]
    c = coords[triangle[3]]

    normal = cross(b - a, c - a)
    offset = -dot(normal, a)

    return Plane(normal, offset)
end

function fundamental_quadratic(planes::Vector{Plane})
    Q = Matrix{Float64}(undef, 4, 4)
    for i in 1:length(planes)
        Q += planes[i].u * planes[i].u'
    end
    return Q
end

function Eh(point::Vector{Float64}, Q::Matrix{Float64})
    return [point;1]' * Q * [point;1]
end

function minEh(Q::Matrix{Float64})
    return Q[1:3, 1:3] \ [0;0;0]
end

# the set of planes spanned by triangles in K incident to at least one vertex in the set Vc
function getIncidentTriangles(K::SimplicialComplex2D, Vc::Set{Int})
    incident_triangles = Set{Triangle}()
    for vertex in Vc
        for edge in K._vertex_to_edges[vertex]
            for triangle in K._edge_to_triangles[edge]
                push!(incident_triangles, triangle)
            end
        end
    end
    return incident_triangles
end

function errorOfEdgeInContractedComplex(K::ContractedSimplicialComplex2D, edge::Edge)
    Q = K._contracted_vertex_Q[edge[1]] + K._contracted_vertex_Q[edge[2]] - K._contracted_edge_Q[edge]
    min_point = minEh(Q)
    return Eh(min_point, Q)
end