using Random
using LinearAlgebra
using Statistics

function triangleEdgesNotSorted(triangle::Triangle)
    return [(triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1])]
end

function rayIntersectsTriangle(ray_origin::Vector{Float32}, ray_direction::Vector{Float32}, v0::Vector{Float32}, v1::Vector{Float32}, v2::Vector{Float32})
    # written by my good friend from india pretending to be an LLM
    EPSILON = 0.000001
    edge1 = v1 - v0
    edge2 = v2 - v0
    h = cross(ray_direction, edge2)
    a = dot(edge1, h)
    if a > -EPSILON && a < EPSILON
        return false
    end
    f = 1.0 / a
    s = ray_origin - v0
    u = f * dot(s, h)
    if u < 0.0 || u > 1.0
        return false
    end
    q = cross(s, edge1)
    v = f * dot(ray_direction, q)
    if v < 0.0 || u + v > 1.0
        return false
    end
    t = f * dot(edge2, q)
    return t > EPSILON
end

function orientedTriangles(sc::SimplicialComplex2D)
    centroid = mean(sc.coords, dims=1)[1]
    all_trig = triangles(sc) # check sorted
    not_oriented_queue = [first(all_trig)]
    oriented_triangles = Dict()

    while length(not_oriented_queue) > 0
        curr_trig = pop!(not_oriented_queue)

        # If triangle is already oriented, skip
        if haskey(oriented_triangles, curr_trig)
            continue
        end

        # Check neighbors and orient the triangles
        for edge in triangleEdges(curr_trig)
            for neigh in sc._edge_to_triangles[edge]
                if neigh != curr_trig && !haskey(oriented_triangles, neigh)
                    push!(not_oriented_queue, neigh)
                end
            end
        end

        # If it's the first triangle, orient it so the normal points away from the centroid
        if length(oriented_triangles) == 0
            oriented_triangles[curr_trig] = curr_trig

            normal_vector = cross((sc.coords[curr_trig[2], :]-sc.coords[curr_trig[1], :])[1], (sc.coords[curr_trig[3], :]-sc.coords[curr_trig[1], :])[1])
            trig_centroid = mean([sc.coords[curr_trig[1], :], sc.coords[curr_trig[2], :], sc.coords[curr_trig[3], :]], dims=1)[1][1]

            # cast a ray from normal vector, count intersections with triangles
            intersections = 0
            for other_trig in all_trig
                if other_trig != curr_trig
                    if rayIntersectsTriangle(trig_centroid, normal_vector, sc.coords[other_trig[1]], sc.coords[other_trig[2]], sc.coords[other_trig[3]])
                        intersections += 1
                    end
                end
            end

            # flip orientation if normal is pointing inwards
            if intersections % 2 == 1
                oriented_triangles[curr_trig] = (curr_trig[1], curr_trig[3], curr_trig[2])
            end
        else
            # Check the orientation of the neighboring triangles
            for edge in triangleEdgesNotSorted(curr_trig)
                for searching_trig in sc._edge_to_triangles[sortEdge(edge)]
                    if searching_trig != curr_trig
                        if haskey(oriented_triangles, searching_trig)
                            # Ensure that the triangle orientation is consistent with its neighbor
                            neighbor_orientation = oriented_triangles[searching_trig]
                            if edge in triangleEdgesNotSorted(neighbor_orientation)
                                oriented_triangles[curr_trig] = (curr_trig[3], curr_trig[2], curr_trig[1])
                            else
                                oriented_triangles[curr_trig] = curr_trig
                            end
                            break
                        end
                    end
                end
            end
        end
    end

    return collect(values(oriented_triangles))
end

function scx2mesh(sc::SimplicialComplex2D)
    vertices = [Point3f(coord...) for coord in sc.coords]
    faces = [TriangleFace(Int.(triangle)...) for triangle in orientedTriangles(sc)]
    return GeometryBasics.normal_mesh(vertices, faces)
end

function visualize(self::ContractedSimplicialComplex2D)
    visualize(self, GeometryBasics.Point{3,Float32}[])
end

function visualize(self::ContractedSimplicialComplex2D, highlights::Array{GeometryBasics.Point{3,Float32}})
    fig = Figure(size=(800, 600))
    ax1 = Axis3(fig[1, 1], aspect=:data)
    ax2 = Axis3(fig[1, 2], aspect=:data)

    orig_mesh = scx2mesh(self.original)
    contracted_mesh = scx2mesh(self.contracted)

    # translation::Vector{Float32}=Vector{Float32}([0.1, 0.0, 0.0])
    # translated_vertices = [v + translation for v in contracted_mesh.position]
    # translated_mesh = GeometryBasics.Mesh(translated_vertices, GeometryBasics.faces(contracted_mesh))

    mesh!(ax1, orig_mesh, color=:red)
    mesh!(ax2, contracted_mesh, color=:blue)

    # hightlight points
    if length(highlights) > 0
        scatter!(ax1, highlights, color=:green)
        scatter!(ax2, highlights, color=:green)
    end

    wireframe!(ax1, orig_mesh, color=:black)
    wireframe!(ax2, contracted_mesh, color=:black)

    display(fig)
end

function visualize(self::SimplicialComplex2D)
    visualize(self, "")
end
function visualize(self::SimplicialComplex2D, highlights::Array{GeometryBasics.Point{3,Float32}}, save_file::String)
    fig = Figure(size=(800, 600))
    ax1 = Axis3(fig[1, 1], aspect=:data)

    orig_mesh = scx2mesh(self)

    mesh!(ax1, orig_mesh, color=:lime)

    # wireframe!(ax1, orig_mesh, color=:black, linewidth=0.5)

    # hightlight points
    if length(highlights) > 0
        scatter!(ax1, highlights, color=:red)
    end

    if save_file != ""
        savefig(fig, save_file)
    else
        display(fig)
    end
end

function visualizeNormals(self::SimplicialComplex2D)
    fig = Figure(size=(800, 600))
    ax1 = Axis3(fig[1, 1], aspect=:data)

    orig_mesh = scx2mesh(self)

    mesh!(ax1, orig_mesh, color=:lime)
    display(fig)
end
