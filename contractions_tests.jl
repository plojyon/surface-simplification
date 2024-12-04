include("contractions.jl")

coords = [Vector([1.0,2.0,3.0]) for i in 1:9]
torus = [(2,1,4), (4,1,3), (5,4,7), (6,5,8), (7,4,6), (7,1,8), (8,2,9), (7,9,1), (1,2,8), (4,5,2), (7,8,5), (3,6,4), (6,9,7), (1,9,3), (8,9,6), (3,2,5), (5,6,3), (2,3,9)]

simplicial_complex = SimplicialComplex2D(torus, coords)
println(simplicial_complex)
