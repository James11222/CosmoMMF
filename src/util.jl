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
    kx = fftfreq(dims[1], sample_rate[1])
    ky = fftfreq(dims[2], sample_rate[2])
    kz = fftfreq(dims[3], sample_rate[3])
    return kx, ky, kz
end

# default T is Float64
wavevectors3D(dims, box_size=(2π, 2π, 2π)) = wavevectors3D(Float64, dims, box_size)
