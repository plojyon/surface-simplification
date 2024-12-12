println("Julia running. Loading packages")

println("DataStructures")
using DataStructures
println("FileIO")
using FileIO
println("Meshes")
using Meshes
println("MeshIO")
using MeshIO
println("GeometryBasics")
using GeometryBasics
println("ProgressBars")
using ProgressBars
println("Serialization")
using Serialization
println("JSON3")
using JSON3
println("Serialization")
using Serialization
println("JLD2")
using JLD2
println("Makie")
using Makie
println("GLMakie")
using GLMakie

println("Including files")
include("contractions.jl")
include("visualisation.jl")
include("preprocessing.jl")

function cachesc(sc, path)
    JSON3.write("$path.json", sc)
    serialize("$path.dat", sc)
    save_object("$path.jld2", sc)
end

function mymain()
    println("Loading bunny")
    #bunidata = loadsc("bunny/reconstruction/bun_zipper.ply")
    bunidata = deserialize("bunidata.dat")

    # rotate by 90 degrees around the x-axis
    for vertex in bunidata.coords
        vertex[2], vertex[3] = vertex[2] * cos(π / 2) - vertex[3] * sin(π / 2), vertex[2] * sin(π / 2) + vertex[3] * cos(π / 2)
    end

    buni = initialContractedSimplicialComplex2D(bunidata)
    fillholes!(buni.contracted)
    calculateFundamentalQuadratics!(buni)

    println("Constructing initial priority queue")
    pq = PriorityQueue()
    for (edge, triangle) in ProgressBar(buni.contracted._edge_to_triangles)
        pq[edge] = error(buni, edge)
    end

    # for _ in ProgressBar(1:34000)
    #     myedge = dequeue!(pq)
    #     while !issafe(buni.contracted, myedge)
    #         myedge = dequeue!(pq)
    #     end

    #     new_vertex = contract!(buni, myedge)
    #     for edge in buni.contracted._vertex_to_edges[new_vertex]
    #         e = error(buni, edge)
    #         pq[edge] = e
    #     end
    # end
    # pq[first(first(pq))]



    fig = Figure(size=(1920, 1080))
    ax = Axis3(fig[1, 1], aspect=:data)

    ax.viewmode = :fit

    start_angle = π / 4
    antispeed = 100

    lastlog = -1000000000
    img_count = 0
    error_hist = []
    record(fig, "bunny_prod.mp4") do io
        tqdm = ProgressBar(1:340000)
        for contraction_count in tqdm
            # if pq[first(first(pq))] > 0.8
            #     break
            # end

            push!(error_hist, pq[first(first(pq))])

            global myedge
            try
                global myedge
                myedge = dequeue!(pq)
                while !issafe(buni.contracted, myedge)
                    myedge = dequeue!(pq)
                end
            catch e
                break
            end

            new_vertex = contract!(buni, myedge)
            for edge in buni.contracted._vertex_to_edges[new_vertex]
                pq[edge] = error(buni, edge)
            end

            thislog = log(length(buni.contracted.vertices))
            if abs(thislog - lastlog) > 0.01
                print('\u0007')
                d_mesh = scx2mesh(buni.contracted)
                meshplt = mesh!(ax, d_mesh, color=:lime)
                wireframeplt = wireframe!(ax, d_mesh, color=:black, linewidth=0.5)
                ax.azimuth[] = start_angle + 2pi * img_count / antispeed
                recordframe!(io)  # record a new frame

                delete!(ax, meshplt)
                delete!(ax, wireframeplt)

                img_count += 1
                lastlog = thislog
            end

            set_description(tqdm, "Frame count: $img_count, thislog: $thislog")
        end
    end

    # plot the error history
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, 1:length(error_hist), error_hist, color=:blue)
    display(fig)
end

mymain()