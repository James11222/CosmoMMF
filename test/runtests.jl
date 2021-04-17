using CosmoMMF
using Test

@testset "CosmoMMF.jl" begin
    # Write your tests here.

function wall(n::Int64)
    """
    Builds a wall out of a nxnxn array. All values in the wall
    are = 1. Other points = 0.
    """
    array = zeros((n,n,n))
    array[:,:,:] .= 0.1
    index = n ÷ 2
    array[index,:,:] .= 1
    return array
end

function cylinder(n::Int64,r::Float64)
    """
    Returns a nxnxn array with a cylinder in the center with
    radius r. All points within the radius are = 1, the rest
    are 0.
    """
    array  = zeros(Float64,(n,n,n))
    array[:,:,:] .= 0.1
    for x::Int64 in 1:n
        for y::Int64 in 1:n
            for z::Int64 in 1:n
                if ((x - (n)/2)^2 + (y - (n)/2)^2) < r^2
                    array[x,y,z] = 1
                end
            end
        end
    end
    return array
end


function sphere(n::Int64,r::Float64)
    """
    Returns a nxnxn array with a sphere in the center with
    radius r. All points within the radius are = 1, the rest
    are 0.
    """
    array = zeros(Float64,(n,n,n))
    array[:,:,:] .= 0.1
    for x::Int64 in 1:n
        for y::Int64 in 1:n
            for z::Int64 in 1:n
                if ((x - (n)/2)^2 + (y - (n)/2)^2 + (z - (n)/2)^2) < r^2
                    array[x,y,z] = 1
                end
            end
        end
    end
    return array
end




end
