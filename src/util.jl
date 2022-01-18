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
    print("This test function works! (totally)")
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

function save_hessian_component(output_directory, Rs, id_string, hessian_component::AbstractArray{T,3}) where T
    """
    A helper function to save a hessian component.
    
    Args:
        output_directory [string]: where should the hessian component be saved to?
        id_string [string]: component number in matrix 
                            -> [[1 2 3],
                                [  4 5],
                                [    6]]
        subvolume_divisions [Int]: number of subvolumes to divide original hessian into on each side. If = 2, then we make 2^3 subvolumes. (Default = 2)
        hessian_component [AbstractArray{T,3}]: one of the components
    """
    dims = size(hessian_component)
    
    div_x, div_y, div_z = dims .÷ 2
    
    #chop hessian into 8ths
    hessian_component_000 = hessian_component[1:div_x, 1:div_y, 1:div_z]
    hessian_component_001 = hessian_component[1:div_x, 1:div_y, div_z+1:dims[3]]
    hessian_component_010 = hessian_component[1:div_x, div_y+1:dims[2], 1:div_z]
    hessian_component_100 = hessian_component[div_x+1:dims[1], 1:div_y, 1:div_z]
    hessian_component_011 = hessian_component[1:div_x, div_y+1:dims[2], div_z+1:dims[3]]
    hessian_component_110 = hessian_component[div_x+1:dims[1], div_y+1:dims[2], 1:div_z]
    hessian_component_101 = hessian_component[div_x+1:dims[1], 1:div_y, div_z+1:dims[3]]
    hessian_component_111 = hessian_component[div_x+1:dims[1], div_y+1:dims[2], div_z+1:dims[3]]
    
    @save output_directory * "hessian_component_000" * "_" * string(Rs) * id_string * ".jld2" hessian_component_000
    @save output_directory * "hessian_component_001" * "_" * string(Rs) * id_string * ".jld2" hessian_component_001
    @save output_directory * "hessian_component_010" * "_" * string(Rs) * id_string * ".jld2" hessian_component_010
    @save output_directory * "hessian_component_100" * "_" * string(Rs) * id_string * ".jld2" hessian_component_100
    @save output_directory * "hessian_component_011" * "_" * string(Rs) * id_string * ".jld2" hessian_component_011
    @save output_directory * "hessian_component_110" * "_" * string(Rs) * id_string * ".jld2" hessian_component_110
    @save output_directory * "hessian_component_101" * "_" * string(Rs) * id_string * ".jld2" hessian_component_101
    @save output_directory * "hessian_component_111" * "_" * string(Rs) * id_string * ".jld2" hessian_component_111
    
end

function compile_sub_hessians(f_Rn::AbstractArray{T,3}, Rs, subvolume, input_directory) where T
    
    id_strings = ["_1", "_2", "_3", "_4", "_5", "_6"]
    dims = size(f_Rn) .÷ 2
    hessian_sub = zeros(Complex{T}, (dims[1], dims[2], dims[3], 6))
    
    for (i, id_string) in enumerate(id_strings)
        if subvolume == "000"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_000
        hessian_sub[:,:,:,i] .= hessian_component_000
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "001"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_001
        hessian_sub[:,:,:,i] .= hessian_component_001
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "010"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_010
        hessian_sub[:,:,:,i] .= hessian_component_010
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "100"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_100
        hessian_sub[:,:,:,i] .= hessian_component_100
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "011"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_011
        hessian_sub[:,:,:,i] .= hessian_component_011
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "110"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_110
        hessian_sub[:,:,:,i] .= hessian_component_110
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "101"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_101
        hessian_sub[:,:,:,i] .= hessian_component_101
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        elseif subvolume == "111"
            @load input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" hessian_component_111
        hessian_sub[:,:,:,i] .= hessian_component_111
            
        #delete old slice we don't need anymore
        file = input_directory * "hessian_component_" * subvolume * "_" * string(Rs) * id_string * ".jld2" 
        run(`rm $file`)
            
        else
            print("Give a correct subvolume string. Example: '000' ")
        end
         
    end
    
    hessian_sub = real(hessian_sub)
    @save input_directory * "hessian_component_" * subvolume * "_full" * ".jld2" hessian_sub   

end
    

function combine_max_sigs(output_directory, subvolumes, field::AbstractArray{T,3}) where T
    dims = size(field)
    div_x, div_y, div_z = size(field) .÷ 2
    sigmax = zeros(T,dims[1],dims[2],dims[3],3)
    
    for subvolume in subvolumes
        @load output_directory * "max_sigs_"*subvolume*".jld2" sigmax_sub
        
        if subvolume == "000"
            sigmax[1:div_x, 1:div_y, 1:div_z, :] .= sigmax_sub
            
        elseif subvolume == "001"
            sigmax[1:div_x, 1:div_y, div_z+1:dims[3], :] .= sigmax_sub
            
        elseif subvolume == "010"
            sigmax[1:div_x, div_y+1:dims[2], 1:div_z, :] .= sigmax_sub
            
        elseif subvolume == "100"
            sigmax[div_x+1:dims[1], 1:div_y, 1:div_z, :] .= sigmax_sub
            
        elseif subvolume == "011"
            sigmax[1:div_x, div_y+1:dims[2], div_z+1:dims[3], :] .= sigmax_sub
            
        elseif subvolume == "110"
            sigmax[div_x+1:dims[1], div_y+1:dims[2], 1:div_z, :] .= sigmax_sub
            
        elseif subvolume == "101"
            sigmax[div_x+1:dims[1], 1:div_y, div_z+1:dims[3], :] .= sigmax_sub
            
        elseif subvolume == "111"
            sigmax[div_x+1:dims[1], div_y+1:dims[2], div_z+1:dims[3], :] .= sigmax_sub
            
        else
            print("Incorrect subvolume string.")
            
        end
            
    end
    
    return sigmax
     
end