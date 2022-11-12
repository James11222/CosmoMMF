using CosmoMMF
using Test

@testset "CosmoMMF.jl" begin
    include("hessian_tests.jl")
    include("filter_tests.jl")
    include("tagging_tests.jl")
end
