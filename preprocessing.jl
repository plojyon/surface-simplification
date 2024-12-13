using NearestNeighbors

function getborder(scx::SimplicialComplex2D)
    border = Set{Edge}()
    for edge in edges(scx)
        if length(scx._edge_to_triangles[edge]) != 2
            push!(border, edge)
        end
    end
    border
end

function getholes(scx::SimplicialComplex2D)
    println("Finding holes")
    holes = Set{Vector{Edge}}()
    border = getborder(scx)
    while length(border) > 0
        println("Holes left: $(length(border))")
        edge = first(border)

        hole = Vector{Edge}([edge]) # edges comprising the hole
        start = edge[1]
        finiš = edge[2]

        # walk around the hole
        while start != finiš
            edges = Set{Edge}(filter(e -> finiš in e, border))
            setdiff!(edges, Set([edge]))
            edge = pop!(edges)
            push!(hole, edge)
            finiš = edge[1] == finiš ? edge[2] : edge[1]
        end
        push!(holes, hole)
        setdiff!(border, hole)
    end
    holes
end

function fillholes!(scx::SimplicialComplex2D)
    holes = getholes(scx)
    println("Filling $(length(holes)) holes")
    for hole in ProgressBar(holes)
        vertices = Vector{Int}()
        for edge in hole
            push!(vertices, edge[1])
            push!(vertices, edge[2])
        end

        # get average coordinate of all vertices
        new_vertex = sum([scx.coords[v] for v in vertices]) / length(vertices)
        push!(scx.coords, new_vertex)
        new_vertex_index = length(scx.coords)
        push!(scx.vertices, new_vertex_index)

        for vertex in vertices
            edge = sortEdge((vertex, new_vertex_index))
            pushToSetDict!(scx._vertex_to_edges, vertex, edge)
            pushToSetDict!(scx._vertex_to_edges, new_vertex_index, edge)
        end
        for edge in hole
            triangle = sortTriangle((new_vertex_index, edge...))
            pushToSetDict!(scx._edge_to_triangles, sortEdge(edge), triangle)
            pushToSetDict!(scx._edge_to_triangles, sortEdge((new_vertex_index, edge[1])), triangle)
            pushToSetDict!(scx._edge_to_triangles, sortEdge((new_vertex_index, edge[2])), triangle)
        end
    end
end

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

    treecoords = Vector{Vector{Float32}}()
    for triangle in mesh
        for vertex in triangle
            push!(treecoords, collect(vertex))
        end
    end
    kdtree = KDTree(hcat(treecoords...); leafsize=25)

    for triangle in ProgressBar(mesh)
        for vertex in triangle
            arr_vertex = collect(vertex)

            # find close neighbours
            idxs = inrange(kdtree, arr_vertex, 0.001)
            neighs = Set{Vector{Float32}}([treecoords[i] for i in idxs])
            intersect!(neighs, keys(anticoords))

            if length(neighs) == 0
                push!(coords, arr_vertex)
                anticoords[arr_vertex] = length(coords)
            else
                anticoords[arr_vertex] = anticoords[first(neighs)]
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
