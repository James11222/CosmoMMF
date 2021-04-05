

"""
    kspace_gaussian_filter(R_S, kv)

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
    for ix = 1:nx, iy = 1:ny, iz = 1:nz
        filter_k[ix, iy, iz] = exp(
            -(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) * R_S^2 / 2)
    end
    return filter_k
end


"""
    kspace_gaussian_filter(f, R_S, kv)

Computes the Fourier representation of a Gaussian filter with width ``R_S``.
You would multiply this with the FFT'd field and then IFFT to apply the filter.

# Arguments:
- `f`: the 3D array to smooth
- `R_S::T`: filter radius
- `kv`: wavevector, tuple of (kx, ky, kz)

# Returns:
- `Array{T,3}`: fourier representation of gaussian filter
"""
function smooth_loggauss(f::AbstractArray{T,3}, R_S::T, kv) where T
    GF = kspace_gaussian_filter(R_S, kv)  # get filter in Fourier space
    f_Rn = 10 .^(
        real(ifft(GF .* fft(log10.(f)))))
    f_Rn = f_Rn .* (sum(f)/sum(f_Rn))  # normalize
    return f_Rn
end


# THESE FUNCTIONS ONLY WORK FOR CUBES. From arxiv:1209.2043

# Computes Hessian of log-gaussian-smoothed density field
function hessian_NEXUSPLUS(f::AbstractArray{T,3}, R_S, kv) where T
    dims = size(f)
    hessian = zeros(Complex{T}, (dims[1], dims[2], dims[3], 6))
    f_Rn = smooth_loggauss(f, R_S, kv)
    f_Rn_hat = fft(f_Rn)

    hescol = 1

    for i in 1:3
        for j in 1:3
            if i <= j # i.e. (1,1), (1,2), (1,3), (2,2), (2,3), (3,3)
                hessian[:,:,:,hescol] = ifft(
                    - kv[i] .* kv[j] .* R_S^2 .* f_Rn_hat )
                hescol += 1
            end
        end
    end
    real(R_S^2 .* hessian)
end


# Computes Hessian of gaussian smoothed density field
function hessian_NEXUS(f::AbstractArray{T,3}, R_S, kv) where T
    dims = size(f)
    hessian = zeros(Complex{T}, (dims[1], dims[2], dims[3], 6))

    f_hat = fft(f)
    filter_times_fhat = exp.( -(kv[1].^2 + kv[2].^2 + kv[3].^2) * R_S^2 / 2  ) .* f_hat

    hescol = 1
    for i in 1:3
        for j in 1:3
            if i <= j # i.e. (1,1), (1,2), (1,3), (2,2), (2,3), (3,3)
                hessian[:,:,:,hescol] = ifft(
                    - kv[i] .* kv[j] .* R_S^2 .* filter_times_fhat
                ) # eq 4 of arxiv:1209.2043
                hescol += 1
            end
        end
    end
    real(R_S^2 * hessian)
end


function signatures_from_hessian(hes::AbstractArray{T,4}) where T
    hsize = size(hes)[1:3]
    cats = zeros(T, (hsize[1], hsize[2], hsize[3], 3) )
    Threads.@threads for i in 1:hsize[1]
        B0 = zeros(T, 3,3)
        eigs = zeros(T, 3)
        hes_slice = zeros(T, 3,3)
        for j in 1:hsize[2]
            for k in 1:hsize[3]

                # StaticArray (SA) here for speed, we know the dimensions are 3x3
                hes_slice = SA[
                    hes[i,j,k,1]  hes[i,j,k,2]  hes[i,j,k,3];
                    hes[i,j,k,2]  hes[i,j,k,4]  hes[i,j,k,5];
                    hes[i,j,k,3]  hes[i,j,k,5]  hes[i,j,k,6]
                ]

                eigs[:] .= sort(real(eigvals(hes_slice)), rev=false)

                cats[i,j,k, 1] = (abs(eigs[3] / eigs[1]) *
                    abs(eigs[3])) * (θ(-eigs[1]) *
                    θ(-eigs[2]) * θ(-eigs[3]))
                cats[i,j,k, 2] = (abs(eigs[2] / eigs[1]) *
                    xθ(1-abs(eigs[3] / eigs[1]))) *
                    (abs(eigs[2]) * θ(-eigs[1]) * θ(-eigs[2]) )
                cats[i,j,k, 3] = (xθ(1-abs(eigs[2] / eigs[1])) *
                    xθ(1-abs(eigs[3] / eigs[1]))) *
                    (abs(eigs[1]) * θ(-eigs[1]))
            end
        end
        GC.gc()
    end
    return cats  # category shape strength ℐ (arxiv:1209.2043 eq 6)
end
