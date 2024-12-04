include("contractions.jl")

function print_struct_properties(obj::T) where T
    println("Properties of $(typeof(obj)):")
    for field in fieldnames(T)
        value = getfield(obj, field)
        println("  $field: $value")
    end
end

coords = [Vector([1.0,2.0,3.0]) for i in 1:9]
torus = [(2,1,4), (4,1,3), (5,4,7), (6,5,8), (7,4,6), (7,1,8), (8,2,9), (7,9,1), (1,2,8), (4,5,2), (7,8,5), (3,6,4), (6,9,7), (1,9,3), (8,9,6), (3,2,5), (5,6,3), (2,3,9)]

simplicial_complex = SimplicialComplex2D(torus, coords)
println("Simplicial complex:")
print_struct_properties(simplicial_complex)

contracted = initialContractedSimplicialComplex2D(simplicial_complex)
println("Initial contracted simplicial complex:")
print_struct_properties(contracted)

println(errorOfEdgeInContractedComplex(contracted, (1,2)))