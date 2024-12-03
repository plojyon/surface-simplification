using DataStructures
using FileIO
using Meshes
using MeshIO
using Makie
using GLMakie

queue = BinaryMinHeap{Int}()


# Pseduocode
# while not isEmpty do ab = MinExtract;
#     if isSafe(ab) then contract ab endif
# endwhile.

function visualize(mesh)
    vertices = [point for triangle in mesh for point in triangle.points]
    vertex_coords = Point3.(vertices)

    faces = Int32[]
    for i in 1:3:length(vertex_coords)
        append!(faces, [i, i + 1, i + 2])
    end
    face_matrix = reshape(faces, :, 3)

    scene = Makie.mesh(vertex_coords, face_matrix, color=:blue)
    display(scene)
end

# buni = load("objects/bunny/reconstruction/bun_zipper.ply")
# visualize(buni)

function isSafe(a, b)
    return true
end

# errorfunction = (a, b) -> a < b  # Min-heap
# ordering = FunctorOrdering(errorfunction)
