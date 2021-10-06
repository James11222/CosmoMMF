
"""
Computes the Fourier representation of a Gaussian filter with width ``R_S``.
You would multiply this with the FFT'd field and then IFFT to apply the filter.

# Arguments:
- `R_S::T`: filter radius
- `kv`: wavevector, tuple of (kx, ky, kz)

# Returns:
- `Array{T,3}`: fourier representation of gaussian filter
"""
function kspace_gaussian_filter(R_S::T, kv) where T
    kx, ky, kz = kv
    nx, ny, nz = length(kx), length(ky), length(kz)
    filter_k = zeros(T, (nx, ny, nz))
    Threads.@threads for ix = 1:nx
        for iy = 1:ny, iz = 1:nz
            filter_k[ix, iy, iz] = exp(
                -(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) * R_S^2 / 2)
        end
    end
    return filter_k
end


"""
Apply a gaussian filter to a density field with smoothing radius `R_S`
making use of wavevectors `kv`

# Arguments:
- `f::AbstractArray{T,3}`: density field
- `R_S::T`: filter radius
- `kv`: wavevector, tuple of (kx, ky, kz)

# Returns:
- Array{T,3}: smoothed density field
"""
function smooth_gauss(f::AbstractArray{T,3}, R_S::T, kv) where T

    GF = kspace_gaussian_filter(R_S, kv)  # get filter in Fourier space
    f_Rn = (real(ifft(GF .* fft(f))))
    f_Rn = f_Rn .* (sum(f)/sum(f_Rn))  # normalize
    return f_Rn
end


"""
Apply a gaussian filter to a density field with smoothing radius `R_S`
making use of wavevectors `kv` using log smoothing

# Arguments:
- `f::AbstractArray{T,3}`: density field
- `R_S::T`: filter radius
- `kv`: wavevector, tuple of (kx, ky, kz)

# Returns:
- Array{T,3}: log-smoothed density field
"""
function smooth_loggauss(f::AbstractArray{T,3}, R_S::T, kv) where {AA, T}
    # get filter in Fourier space
    GF = kspace_gaussian_filter(R_S, kv)

    # # simple version of what follows
    # f_Rn = 10 .^(real(ifft(GF .* fft(log10.(f)))))

    buffer = similar(f, Complex{T})   # make a complex version of the array
    @strided buffer .= log10.(f)      # threaded log10 the entire array
    fft!(buffer)                      # FFT
    @strided buffer .= GF .* buffer   # convolve
    ifft!(buffer)                     # IFFT
    f_result = real(buffer)           # get the real part
    @strided f_result .= (10 .^ f_result)  # threaded exp
    norm = sum(f)/sum(f_result)       # compute normalization
    @strided f_result .= f_result .* norm  # normalize
    return f_result
end



# THESE FUNCTIONS ONLY WORK FOR CUBES. From arxiv:1209.2043

# Computes Hessian of log-gaussian-smoothed density field
# f_Rn is the smoothed field
function hessian_from_smoothed(f_Rn::AbstractArray{T,3}, R_S, kv) where T
    dims = size(f_Rn)
    hessian = zeros(Complex{T}, (dims[1], dims[2], dims[3], 6))
    f_Rn_hat = fft(f_Rn)

    kx = kv[1]
    ky = kv[2]
    kz = kv[3]

    Threads.@threads for x in 1:length(kv[1])
        for y in 1:length(kv[2])
            for z in 1:length(kv[3])
                # (1,1)
                hessian[x,y,z,1] =
                    - kx[x] * kx[x] * R_S^2 * f_Rn_hat[x,y,z]
                # (1,2)
                hessian[x,y,z,2] =
                    - kx[x] * ky[y] * R_S^2 * f_Rn_hat[x,y,z]
                # (1,3)
                hessian[x,y,z,3] =
                    - kx[x] * kz[z] * R_S^2 * f_Rn_hat[x,y,z]
                # (2,2)
                hessian[x,y,z,4] =
                    - ky[y] * ky[y] * R_S^2 * f_Rn_hat[x,y,z]
                # (2,3)
                hessian[x,y,z,5] =
                    - ky[y] * kz[z] * R_S^2 * f_Rn_hat[x,y,z]
                # (3,3)
                hessian[x,y,z,6] =
                    - kz[z] * kz[z] * R_S^2 * f_Rn_hat[x,y,z]
            end
        end
    end
    hessian[:,:,:,1] .= ifft(hessian[:,:,:,1])
    hessian[:,:,:,2] .= ifft(hessian[:,:,:,2])
    hessian[:,:,:,3] .= ifft(hessian[:,:,:,3])
    hessian[:,:,:,4] .= ifft(hessian[:,:,:,4])
    hessian[:,:,:,5] .= ifft(hessian[:,:,:,5])
    hessian[:,:,:,6] .= ifft(hessian[:,:,:,6])
    real(R_S^2 .* hessian)
end


function signatures_from_hessian(hes::AbstractArray{T,4}) where T
    hsize = size(hes)[1:3]
    sigs = zeros(T, (hsize[1], hsize[2], hsize[3], 3) )
    Threads.@threads for i in 1:hsize[1]
        for j in 1:hsize[2]
            for k in 1:hsize[3]

                # StaticArray (SA) here for speed, we know the dimensions are 3x3
                #the hessian matrix is symmetric
                hes_slice = SA[
                    hes[i,j,k,1]  hes[i,j,k,2]  hes[i,j,k,3];
                    hes[i,j,k,2]  hes[i,j,k,4]  hes[i,j,k,5];
                    hes[i,j,k,3]  hes[i,j,k,5]  hes[i,j,k,6]
                ]

                e1, e2, e3 = sort(real(eigvals(hes_slice)), rev=false)

                sigs[i,j,k,1] = (abs(e3 / e1) *
                    abs(e3)) * (θ(-e1) * θ(-e2) * θ(-e3))
                sigs[i,j,k,2] = (abs(e2 / e1) *
                    xθ(1-abs(e3 / e1))) * (abs(e2) * θ(-e1) * θ(-e2) )
                sigs[i,j,k,3] = (xθ(1-abs(e2 / e1)) *
                    xθ(1-abs(e3 / e1))) * (abs(e1) * θ(-e1))
            end
        end
    end

    # category shape strength 𝒮 (arxiv:1209.2043 eq 7)
    return sigs
end


"""
In this function we start at minimum length scale to calculate signatures
then iterate to the maximum scale of the problem, calculating the signatures
for each smoothing scale. We then chose the maximum signature over all scales
at each point to categorize whether that point is a wall, filament, or cluster.
"""
function (Rs, field::AbstractArray{T,3}; alg=:NEXUSPLUS) where T

    if alg ∉ (:NEXUS, :NEXUSPLUS)
        throw(ArgumentError("alg must be either :NEXUS or :NEXUSPLUS"))
    end

    nx,ny,nz = size(field)
    
    #we need to make sure the field has no 0 values
    field = field .+ 0.0001

    #calculate wave vectors for our field
    wave_vecs = wavevectors3D(T, (nx,ny,nz))

    #calculate signatures at each scale Rn, determine max
    sigmax = zeros(T,nx,ny,nz,3)

    for (R_index, R) in enumerate(Rs)

        if alg == :NEXUS
            f_Rn = smooth_gauss(field, R, wave_vecs)
        elseif alg == :NEXUSPLUS
            f_Rn = smooth_loggauss(field, R, wave_vecs)
        end
#         print("succesfully calculated f_Rn \n")

        H_Rn = hessian_from_smoothed(f_Rn, R, wave_vecs)
        
#         print("Succesfully calculated H_Rn \n")

        sigs_Rn = signatures_from_hessian(H_Rn)
        
#         print("Succesfully calculated S_Rn \n")
        
        for ix = 1:nx, iy = 1:ny, iz = 1:nz, sigtype = 1:3
            sigmax[ix, iy, iz, sigtype] = max(
                sigmax[ix, iy, iz, sigtype], sigs_Rn[ix, iy, iz, sigtype])
        end
    end

    return sigmax
end
