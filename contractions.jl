using LinearAlgebra

const Triangle = Tuple{Int,Int,Int}
const Edge = Tuple{Int,Int}

function triangleEdges(triangle::Triangle)
    return [(triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[1], triangle[3])]
end

function sortEdge(edge::Edge)
    return edge[1] < edge[2] ? edge : (edge[2], edge[1])
end

function pushToSetDict!(dict::Dict{K,Set{V}}, key::K, value::V) where {K,V}
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = Set([value])
    end
end

function sortTriangle(triangle::Triangle)
    if triangle[1] < triangle[2] < triangle[3]
        return triangle
    elseif triangle[1] < triangle[3] < triangle[2]
        return (triangle[1], triangle[3], triangle[2])
    elseif triangle[2] < triangle[1] < triangle[3]
        return (triangle[2], triangle[1], triangle[3])
    elseif triangle[2] < triangle[3] < triangle[1]
        return (triangle[2], triangle[3], triangle[1])
    elseif triangle[3] < triangle[1] < triangle[2]
        return (triangle[3], triangle[1], triangle[2])
    else
        return (triangle[3], triangle[2], triangle[1])
    end
end

struct SimplicialComplex2D
    coords::Array{Vector{Float32}}
    vertices::Set{Int}

    _vertex_to_edges::Dict{Int,Set{Edge}}
    _edge_to_triangles::Dict{Edge,Set{Triangle}}

    function SimplicialComplex2D(triangles::Array{Triangle}, coords::Array{Vector{Float32}})
        _vertex_to_edges = Dict{Int,Set{Edge}}()
        _edge_to_triangles = Dict{Edge,Set{Triangle}}()

        for triangle in triangles
            triangle = sortTriangle(triangle)
            for edge in triangleEdges(triangle)
                edge = sortEdge(edge)
                pushToSetDict!(_edge_to_triangles, edge, triangle)

                for vertex in edge
                    pushToSetDict!(_vertex_to_edges, vertex, edge)
                end
            end
        end

        new(coords, keys(_vertex_to_edges), _vertex_to_edges, _edge_to_triangles)
    end
end

function edges(K::SimplicialComplex2D)
    out = Set(union(values(K._vertex_to_edges)...))
end

function triangles(K::SimplicialComplex2D)
    return Set(union(values(K._edge_to_triangles)...))
end

function removeEdge(K::SimplicialComplex2D, edge::Edge)
    delete!(K._vertex_to_edges[edge[1]], edge)
    delete!(K._vertex_to_edges[edge[2]], edge)
end

function removeTriangle(K::SimplicialComplex2D, triangle::Triangle)
    for edge in [(triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1])]
        delete!(K._edge_to_triangles[edge], triangle)
    end
end

# replaces the old_vertex with the new_vertex in the simplicial complex K
# does not handle adding new coords
function replaceVertexInEdgesAndTriangles(K::SimplicialComplex2D, old_vertex::Int, new_vertex::Int)
    delete!(K.vertices, old_vertex)
    push!(K.vertices, new_vertex)

    for edge in K._vertex_to_edges[old_vertex]
        # construct the new edge with the new vertex and remove the old one
        if edge[1] == old_vertex
            new_edge = (new_vertex, edge[2])

            # remove the edge from the other vertex that is not being replaced
            delete!(K._vertex_to_edges[edge[2]], edge)
        else
            new_edge = (edge[1], new_vertex)
            delete!(K._vertex_to_edges[edge[1]], edge)
        end

        new_edge = sortEdge(new_edge)

        # add the new edge to the vertex_to_edges dictionary
        pushToSetDict!(K._vertex_to_edges, new_edge[1], new_edge)
        pushToSetDict!(K._vertex_to_edges, new_edge[2], new_edge)

        for triangle in K._edge_to_triangles[edge]
            new_triangle = sortTriangle(Tuple([vertex == old_vertex ? new_vertex : vertex for vertex in triangle]))

            # triangle is stored in 3 entries edges => triangles
            for edge_of_trig in triangleEdges(triangle)
                # we do not need to delete the triangle from the old edge
                # because the old edge will be removed as a key
                if edge_of_trig != edge
                    delete!(K._edge_to_triangles[edge_of_trig], triangle)
                    pushToSetDict!(K._edge_to_triangles, new_edge, new_triangle)
                end
            end
        end
        delete!(K._edge_to_triangles, edge)
    end

    delete!(K._vertex_to_edges, old_vertex)
end

struct ContractedSimplicialComplex2D
    original::SimplicialComplex2D
    contracted::SimplicialComplex2D

    _contracted_triangle_Q::Dict{Triangle,Matrix{Float32}}
    _contracted_edge_Q::Dict{Edge,Matrix{Float32}}
    _contracted_vertex_Q::Dict{Int,Matrix{Float32}}

    function ContractedSimplicialComplex2D(original, contracted)
        _contracted_triangle_Q = Dict{Triangle,Matrix{Float32}}()
        _contracted_edge_Q = Dict{Edge,Matrix{Float32}}()
        _contracted_vertex_Q = Dict{Int,Matrix{Float32}}()

        new(original, contracted, _contracted_triangle_Q, _contracted_edge_Q, _contracted_vertex_Q)
    end
end

function initialContractedSimplicialComplex2D(K::SimplicialComplex2D)
    contracted = ContractedSimplicialComplex2D(K, deepcopy(K))

    calculateFundamentalQuadratics(contracted)
    return contracted
end

# fills in the _triangle_Q, _edge_Q, and _vertex_Q fields of K
function calculateFundamentalQuadratics(K::ContractedSimplicialComplex2D)
    for triangle in triangles(K.contracted)
        K._contracted_triangle_Q[triangle] = fundamental_quadratic([planeFromTriangle(triangle, K.contracted.coords)])
    end

    for edge in edges(K.contracted)
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
    normal::Vector{Float32}
    offset::Float32
    u::Vector{Float32}

    function Plane(normal::Vector{Float32}, offset::Float32)
        new(normal, offset, [normal; offset])
    end
end

function planeFromTriangle(triangle::Triangle, coords::Array{Vector{Float32}})
    #TODO [CHATGPT] confirm code

    a = coords[triangle[1]]
    b = coords[triangle[2]]
    c = coords[triangle[3]]

    normal = cross(b - a, c - a)
    offset = -dot(normal, a)

    return Plane(normal, offset)
end

function fundamental_quadratic(planes::Vector{Plane})
    Q = Matrix{Float32}(undef, 4, 4)
    for i in 1:length(planes)
        Q += planes[i].u * planes[i].u'
    end
    return Q
end

function Eh(point::Vector{Float32}, Q::Matrix{Float32})
    return [point; 1]' * Q * [point; 1]
end

function minEh(Q::Matrix{Float32})
    return Q[1:3, 1:3] \ [0; 0; 0]
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

function minCForEdgeContraction(K::ContractedSimplicialComplex2D, edge::Edge)
    Q = K._contracted_vertex_Q[edge[1]] + K._contracted_vertex_Q[edge[2]] - K._contracted_edge_Q[edge]
    return minEh(Q)
end

function contract!(K::ContractedSimplicialComplex2D, edge::Edge)
    # 1) remove ab, abx, and aby
    # remove edge
    delete!(K.contracted._vertex_to_edges[edge[1]], edge)
    delete!(K.contracted._vertex_to_edges[edge[2]], edge)

    # remove triangles (each triangle is stored in 3 entries edges => triangles)
    for triangle in K.contracted._edge_to_triangles[edge]
        for edge_trig in triangleEdges(triangle)
            delete!(K.contracted._edge_to_triangles[edge_trig], triangle)
        end
    end
    delete!(K.contracted._edge_to_triangles, edge)

    # 2) substitute c for a and for b wherever they occur in the remaining set of
    # vertices, edges, and triangles (removing resulting duplications making sure L is a set.)
    c_coords = minCForEdgeContraction(K, edge)
    c = length(K.contracted.coords) + 1
    push!(K.contracted.coords, c_coords)

    replaceVertexInEdgesAndTriangles(K.contracted, edge[1], c)
    replaceVertexInEdgesAndTriangles(K.contracted, edge[2], c)
end
