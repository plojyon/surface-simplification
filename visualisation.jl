
function scx2mesh(sc::SimplicialComplex2D)
    vertices = [Point3f(coord...) for coord in sc.coords]
    faces = [TriangleFace(Int.(triangle)...) for triangle in triangles(sc)]
    return GeometryBasics.Mesh(vertices, faces)
end

function visualize(self::SimplicialComplex2D)
    Makie.mesh(tomesh(self))
end

function visualize(self::ContractedSimplicialComplex2D)
    visualize(self, GeometryBasics.Point{3,Float32}[])
end

function visualize(self::ContractedSimplicialComplex2D, highlights::Array{GeometryBasics.Point{3,Float32}})
    fig = Figure(size=(800, 600))
    ax1 = Axis3(fig[1, 1])
    ax2 = Axis3(fig[1, 2])

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

    # wireframe!(ax1, orig_mesh, color=:black)
    # wireframe!(ax2, contracted_mesh, color=:black)

    display(fig)
end
