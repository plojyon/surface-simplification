struct SimplicialComplex2D
    coords::Array{Vector{Float64}}
    vertices::Set{Int}
    edges::Set{Tuple{Int,Int}}
    triangles::Set{Tuple{Int,Int,Int}}

    _vertex_to_edges::Dict{Int,Vector{Int}}
    _edge_to_triangles::Dict{Int,Vector{Int}}

    function SimplicialComplex2D(triangles::Vector{Tuple{Float64,Float64,Float64}})
        coords = Vector{Vector{Float64}}()
        vertices = Set{Int}()
        edges = Set{Tuple{Int,Int}}()
        triangles = Set{Tuple{Int,Int,Int}}()

        _vertex_to_edges = Dict{Int,Vector{Int}}()
        _edge_to_triangles = Dict{Int,Vector{Int}}()

        for triangle in triangles
            for vertex in triangle
                push!(vertices, vertex)
            end
            push!(triangles, triangle)
        end

        for vertex in vertices
            push!(coords, [0.0, 0.0, 0.0])
        end

        for triangle in triangles
            for i in 1:3
                edge = (triangle[i], triangle[i % 3 + 1])
                if edge[1] > edge[2]
                    edge = (edge[2], edge[1])
                end
                push!(edges, edge)
            end
        end

        for edge in edges
            push!(_vertex_to_edges[edge[1]], edge)
            push!(_vertex_to_edges[edge[2]], edge)
        end

        for triangle in triangles
            for i in 1:3
                edge = (triangle[i], triangle[i % 3 + 1])
                if edge[1] > edge[2]
                    edge = (edge[2], edge[1])
                end
                push!(_edge_to_triangles[edge], triangle)
            end
        end

        new(coords, vertices, edges, triangles, _vertex_to_edges, _edge_to_triangles)
    end
end

struct ContractedSimplicialComplex2D
    original::SimplicialComplex2D
    contracted::SimplicialComplex2D
    vertex_contracted_to_original::Dict{Int,Set{Int}}

    _contracted_triangle_Q::Dict{Tuple{Int,Int,Int},Matrix{Float64}} = Dict()
    _contracted_edge_Q::Dict{Tuple{Int,Int},Matrix{Float64}} = Dict()
    _contracted_vertex_Q::Dict{Int,Matrix{Float64}} = Dict()
end

# fills in the _triangle_Q, _edge_Q, and _vertex_Q fields of K
function calculateFundamentalQuadratics(K::ContractedSimplicialComplex2D)
    for triangle in K.contracted.triangles
        K._contracted_triangle_Q[triangle] = fundamental_quadratic([planeFromTriangle(triangle, K.contracted.coords)])
    end

    for edge in K.contracted.edges
        neigh_triangles = K.contracted._edge_to_triangles[edge]
        K._contracted_edge_Q[edge] = sum([K._contracted_edge_Q[triangle] for triangle in neigh_triangles])
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

function planeFromTriangle(triangle::Tuple{Int,Int,Int}, coords::Array{Vector{Float64}})
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

function distance(point::Vector{Float64}, Q::Matrix{Float64})
    return point' * Q * point
end

function minDistancePoint(point::Vector{Float64}, Q::Matrix{Float64})
    sol = Q \ [point; 1.0]
    return sol[1:3]
end

# the set of planes spanned by triangles in K incident to at least one vertex in the set Vc
function getIncidentTriangles(K::SimplicialComplex2D, Vc::Set{Int})
    incident_triangles = Set{Int}()
    for vertex in Vc
        for edge in K._vertex_to_edges[vertex]
            for triangle in K._edge_to_triangles[edge]
                push!(incident_triangles, triangle)
            end
        end
    end
    return incident_triangles
end

function errorOfEdgeInContractedComplex(K::ContractedSimplicialComplex2D, edge::Tuple{Int,Int})
    # Vc is the set of all points that were in K0 and contracted to
    # c where c is the contracted edge
    Vc = K.vertex_contracted_to_original[edge[1]] ∪ K.vertex_contracted_to_original[edge[2]] ∪ Set([edge[1], edge[2]])

    # This is labeled as H in the paper
    incident_triangles = getIncidentTriangles(K.original, Vc)

    Q = getFundamentalQuadratic(K, incident_triangles)
end