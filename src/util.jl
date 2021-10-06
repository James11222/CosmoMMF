# computational utilities

"""
    θ(x)

Heaviside step function. Returns one if the argument is positive, zero otherwise.
"""
function θ(x::T) where T
    if x > zero(T)
        return one(T)
    end
    return zero(T)
end

"""
    xθ(x)

Heaviside step times x. Returns x if the argument is positive, zero otherwise.
"""
function xθ(x::T) where T
    if x > zero(T)
        return x
    end
    return zero(T)
end


"""
    wavevectors3D([T=Float64], dims, box_size=(2π, 2π, 2π))

# Arguments:
- `T` : output element type
- `dims`: array or tuple containing size
- `box_size`: side lengths of box, normalizes the wavenumber

# Returns:
- `Tuple{Vector{T},Vector{T},Vector{T}}`: kx, ky, kz wavevectors

# Examples
```julia-repl
julia> CosmoMMF.wavevectors3D((6,2,2))
([0.0, 1.0, 2.0, -3.0, -2.0, -1.0], [0.0, -1.0], [0.0, -1.0])
```
"""
function wavevectors3D(T::Type{<:Real}, dims, box_size=(2π, 2π, 2π))
    sample_rate = T.(2π .* dims ./ box_size)
    kx = fftfreq(dims[1], sample_rate[1]) .* (2π / dims[1])
    ky = fftfreq(dims[2], sample_rate[2]) .* (2π / dims[2])
    kz = fftfreq(dims[3], sample_rate[3]) .* (2π / dims[3])
    return kx, ky, kz
end

# default T is Float64
wavevectors3D(dims, box_size=(2π, 2π, 2π)) = wavevectors3D(Float64, dims, box_size)


function test_print() 
    print("This test function works!")
end

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


function save_max_sigs(output_directory, id_string, sig_array::AbstractArray{T,4}) where T
    """
    A simple function to save the signatures to an output directory in the form of a jld2 file.
    """
    @save output_directory * id_string * ".jld2" sig_array
end
    
