function kspace_gaussian_filter(R_S::T, kv) where T
    """
    Computes the Fourier representation of a Gaussian filter with width ``R_S``.
    You would multiply this with the FFT'd field and then IFFT to apply the filter.

    # Arguments:
    - `R_S::T`: filter radius
    - `kv`: wavevector, tuple of (kx, ky, kz)

    # Returns:
    - `Array{T,3}`: fourier representation of gaussian filter
    """
    kx, ky, kz = kv
    nx, ny, nz = length(kx), length(ky), length(kz)
    filter_k = zeros(T, (nx, ny, nz))
    for ix = 1:nx, iy = 1:ny, iz = 1:nz
        filter_k[ix, iy, iz] = exp(
            -(kx[ix]^2 + ky[iy]^2 + kz[iz]^2) * R_S^2 / 2)
    end
    return filter_k
end

function smooth_gauss(f::AbstractArray{T,3}, R_S::T, kv) where T
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
    GF = kspace_gaussian_filter(R_S, kv)  # get filter in Fourier space
    f_Rn = (real(ifft(GF .* fft(f))))
    f_Rn = f_Rn .* (sum(f)/sum(f_Rn))  # normalize
    return f_Rn
end

function smooth_loggauss(f::AbstractArray{T,3}, R_S::T, kv) where T
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
    
    kx = kv[1]
    ky = kv[2]
    kz = kv[3]
    
    for x in 1:length(kv[1])
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
    hessian[:,:,:,1] = ifft(hessian[:,:,:,1])
    hessian[:,:,:,2] = ifft(hessian[:,:,:,2])
    hessian[:,:,:,3] = ifft(hessian[:,:,:,3])
    hessian[:,:,:,4] = ifft(hessian[:,:,:,4])
    hessian[:,:,:,5] = ifft(hessian[:,:,:,5])
    hessian[:,:,:,6] = ifft(hessian[:,:,:,6])
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
    sigs = zeros(T, (hsize[1], hsize[2], hsize[3], 3) )
    Threads.@threads for i in 1:hsize[1]
        B0 = zeros(T, 3,3)
        eigs = zeros(T, 3)
        hes_slice = zeros(T, 3,3)
        for j in 1:hsize[2]
            for k in 1:hsize[3]

                # StaticArray (SA) here for speed, we know the dimensions are 3x3
                #the hessian matrix is symmetric
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
                
                sigs[i,j,k,1] = cats[i,j,k, 1] * abs(eigs[3])*θ(-eigs[1])*θ(-eigs[2])*θ(-eigs[3])
                sigs[i,j,k,2] = cats[i,j,k, 2] * abs(eigs[2])*θ(-eigs[1])*θ(-eigs[2])
                sigs[i,j,k,3] = cats[i,j,k, 3] * abs(eigs[1])*θ(-eigs[1])
            end
        end
        # GC.gc()
    end
    #return cats  # category shape strength ℐ (arxiv:1209.2043 eq 6)
    return sigs
end


function maximum_signature(R_0::Float64, field::AbstractArray{T,3}) where T
    """
    In this function we start at minimum length scale to calculate signatures
    then iterate to the maximum scale of the problem, calculating the signatures
    for each smoothing scale. We then chose the maximum signature over all scales
    at each point to categorize whether that point is a wall, filament, or cluster.

    R_0 is typically the grid spacing of the field
    """

    nx,ny,nz = size(field)

    #create an array of scales
    Rs = []
    n::Int64 = 5 #this will need to be changed
    for i in 1:n
        Ri = R_0 * sqrt(2)^i
        Rs = append!(Rs,Ri)
    end

    #calculate wave vectors for our field
    wave_vecs = wavevectors3D((nx,ny,nz))

    #calculate signatures at each scale Rn then determine the max signature
    sigs_ = zeros(T,nx,ny,nz,3,length(Rs))
    sigs_final = zeros(T,nx,ny,nz,3)

    for a in 1:length(Rs)
        log_smooth_field = smooth_loggauss(field, Rs[a], wave_vecs)
        H_Rn = hessian_NEXUSPLUS(log_smooth_field, Rs[a], wave_vecs)
        sigs_Rn = signatures_from_hessian(H_Rn)
        sig_c,sig_f,sig_w = sigs_Rn[:,:,:,1],sigs_Rn[:,:,:,2],sigs_Rn[:,:,:,3]
        sigs_[:,:,:,1,a] .= sig_c
        sigs_[:,:,:,2,a] .= sig_f
        sigs_[:,:,:,3,a] .= sig_w
    end

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                sig_max_per_Rn_cluster = maximum(sigs_[i,j,k,1,:])
                sig_max_per_Rn_filament = maximum(sigs_[i,j,k,2,:])
                sig_max_per_Rn_wall = maximum(sigs_[i,j,k,3,:])
                sigs_final[i,j,k,1] = sig_max_per_Rn_cluster
                sigs_final[i,j,k,2] = sig_max_per_Rn_filament
                sigs_final[i,j,k,3] = sig_max_per_Rn_wall
            end
        end
    end

    return sigs_final
end
