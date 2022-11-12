using Test 
using CosmoMMF
using JLD2

@testset "Test_Hessian" begin
    
    resolution = 64
    Rs = (√2) .^ 4
    wave_vecs = CosmoMMF.wavevectors3D((resolution,resolution,resolution))
    
    #generate some test fields
    sphere_field = CosmoMMF.sphere(resolution, 5.0)
    smooth_sphere_field = CosmoMMF.smooth_gauss(sphere_field, Rs, wave_vecs)
    
    #compute hessians
    H_slow = CosmoMMF.slow_hessian_from_smoothed(smooth_sphere_field, Rs, wave_vecs)
    H_fast = CosmoMMF.fast_hessian_from_smoothed(smooth_sphere_field, Rs, wave_vecs)

    @test H_slow == H_fast
    
end