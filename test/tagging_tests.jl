using Test 
using CosmoMMF
using JLD2


@testset "Test_Tagging" begin
    Rs = (√2) .^ 0:5

    @load "test_data/test_clusbool.jld2" clusbool
    @load "test_data/test_density_field.jld2" den

    #compute maximum structure signatures
    max_sigs_nexusplus = CosmoMMF.maximum_signature(Rs, den, alg=:NEXUSPLUS)

    #tagging
    clusbool, filbool, wallbool, voidbool, S_fil, dM2_fil, S_wall, dM2_wall = CosmoMMF.calc_structure_bools(
    den, max_sigs_nexusplus, true, clusbool)

    #test to make sure the filament curve is well behaved
    fil_max_index = findmax(dM2_fil ./ maximum(dM2_fil))[2]
    fil_test_passed = dM2_fil[fil_max_index - 1] < dM2_fil[fil_max_index] || dM2_fil[fil_max_index + 1] < dM2_fil[fil_max_index]
    @test fil_test_passed

    #test to make sure the wall curve is well behaved
    wall_max_index = findmax(dM2_wall ./ maximum(dM2_wall))[2]
    wall_test_passed = dM2_wall[wall_max_index - 1] < dM2_wall[wall_max_index] || dM2_wall[wall_max_index + 1] < dM2_wall[wall_max_index]
    @test wall_test_passed

    #test the virialization technique
    clusbool, filbool, wallbool, voidbool, S_fil, dM2_fil, S_wall, dM2_wall = CosmoMMF.calc_structure_bools(
    den, max_sigs_nexusplus, true)

    #test to make sure the filament curve is well behaved
    fil_max_index = findmax(dM2_fil ./ maximum(dM2_fil))[2]
    fil_test_passed = dM2_fil[fil_max_index - 1] < dM2_fil[fil_max_index] || dM2_fil[fil_max_index + 1] < dM2_fil[fil_max_index]
    @test fil_test_passed

    #test to make sure the wall curve is well behaved
    wall_max_index = findmax(dM2_wall ./ maximum(dM2_wall))[2]
    wall_test_passed = dM2_wall[wall_max_index - 1] < dM2_wall[wall_max_index] || dM2_wall[wall_max_index + 1] < dM2_wall[wall_max_index]
    @test wall_test_passed

end