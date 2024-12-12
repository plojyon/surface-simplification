using Random
using LinearAlgebra

function triangleEdgesNotSorted(triangle::Triangle)
    return [(triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1])]
end

function orientedTriangles(sc::SimplicialComplex2D)
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

        # If it's the first triangle, we can orient it arbitrarily
        if length(oriented_triangles) == 0
            oriented_triangles[curr_trig] = (curr_trig[3], curr_trig[2], curr_trig[1])
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
function visualize(self::SimplicialComplex2D, save_file::String)
    fig = Figure(size=(800, 600))
    ax1 = Axis3(fig[1, 1], aspect=:data)

    orig_mesh = scx2mesh(self)

    mesh!(ax1, orig_mesh, color=:lime)

    # wireframe!(ax1, orig_mesh, color=:black, linewidth=0.5)

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