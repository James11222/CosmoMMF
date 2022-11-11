using Test 
using CosmoMMF
using JLD2


@testset "Test_Max_Signatures" begin
    resolution = 64

    #generate some test fields
    sphere_field = CosmoMMF.sphere(resolution,5.0)
    cylinder_field = CosmoMMF.cylinder(resolution,5.0)
    wall_field = CosmoMMF.wall(resolution)

    Rs = (√2) .^ 5

    #calculate max signatures
    max_sigs_sphere = CosmoMMF.maximum_signature(Rs, sphere_field)
    max_sigs_cylinder = CosmoMMF.maximum_signature(Rs, cylinder_field)
    max_sigs_wall = CosmoMMF.maximum_signature(Rs, wall_field)

    #we compare the simple implementation to a method that reduces RAM footprint
    output_directory = "temp_output/"

    save_name = "test_sphere"
    max_sigs_sphere_RR = CosmoMMF.reduce_RAM_maximum_signature(Rs, output_directory, save_name, sphere_field)
    @load "temp_output/max_sigs_test_sphere_full_signatures.jld2" sigmax
    max_sigs_sphere_RR = sigmax

    save_name = "test_cylinder"
    max_sigs_cylinder_RR = CosmoMMF.reduce_RAM_maximum_signature(Rs, output_directory, save_name, cylinder_field)
    @load "temp_output/max_sigs_test_cylinder_full_signatures.jld2" sigmax
    max_sigs_cylinder_RR = sigmax

    save_name = "test_wall"
    max_sigs_wall_RR = CosmoMMF.reduce_RAM_maximum_signature(Rs, output_directory, save_name, wall_field)
    @load "temp_output/max_sigs_test_wall_full_signatures.jld2" sigmax
    max_sigs_wall_RR = sigmax

    #we check to see if the two methods are equivalent where the structures were originally located
    bool_filter = sphere_field .== 1
    @test max_sigs_sphere[bool_filter,:] == max_sigs_sphere_RR[bool_filter, :]

    bool_filter = cylinder_field .== 1
    @test max_sigs_cylinder[bool_filter,:] == max_sigs_cylinder_RR[bool_filter, :] 
    
    bool_filter = wall_field .== 1
    @test max_sigs_wall[bool_filter,:] == max_sigs_wall_RR[bool_filter, :]

end