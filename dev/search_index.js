var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = CosmoMMF","category":"page"},{"location":"#CosmoMMF","page":"Home","title":"CosmoMMF","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [CosmoMMF]","category":"page"},{"location":"#CosmoMMF.kspace_gaussian_filter-Union{Tuple{T}, Tuple{T, Any}} where T","page":"Home","title":"CosmoMMF.kspace_gaussian_filter","text":"Computes the Fourier representation of a Gaussian filter with width R_S. You would multiply this with the FFT'd field and then IFFT to apply the filter.\n\nArguments:\n\nR_S::T: filter radius\nkv: wavevector, tuple of (kx, ky, kz)\n\nReturns:\n\nArray{T,3}: fourier representation of gaussian filter\n\n\n\n\n\n","category":"method"},{"location":"#CosmoMMF.maximum_signature-Union{Tuple{T}, Tuple{Any, AbstractArray{T, 3}}} where T","page":"Home","title":"CosmoMMF.maximum_signature","text":"In this function we start at minimum length scale to calculate signatures then iterate to the maximum scale of the problem, calculating the signatures for each smoothing scale. We then chose the maximum signature over all scales at each point to categorize whether that point is a wall, filament, or cluster.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoMMF.reduce_RAM_maximum_signature-Union{Tuple{T}, Tuple{Any, Any, Any, AbstractArray{T, 3}}} where T","page":"Home","title":"CosmoMMF.reduce_RAM_maximum_signature","text":"In this function we start at minimum length scale to calculate signatures then iterate to the maximum scale of the problem, calculating the signatures for each smoothing scale. We then chose the maximum signature over all scales at each point to categorize whether that point is a wall, filament, or cluster.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoMMF.smooth_gauss-Union{Tuple{T}, Tuple{AbstractArray{T, 3}, T, Any}} where T","page":"Home","title":"CosmoMMF.smooth_gauss","text":"Apply a gaussian filter to a density field with smoothing radius R_S making use of wavevectors kv\n\nArguments:\n\nf::AbstractArray{T,3}: density field\nR_S::T: filter radius\nkv: wavevector, tuple of (kx, ky, kz)\n\nReturns:\n\nArray{T,3}: smoothed density field\n\n\n\n\n\n","category":"method"},{"location":"#CosmoMMF.smooth_loggauss-Union{Tuple{T}, Tuple{AA}, Tuple{AbstractArray{T, 3}, T, Any}} where {AA, T}","page":"Home","title":"CosmoMMF.smooth_loggauss","text":"Apply a gaussian filter to a density field with smoothing radius R_S making use of wavevectors kv using log smoothing\n\nArguments:\n\nf::AbstractArray{T,3}: density field\nR_S::T: filter radius\nkv: wavevector, tuple of (kx, ky, kz)\n\nReturns:\n\nArray{T,3}: log-smoothed density field\n\n\n\n\n\n","category":"method"},{"location":"#CosmoMMF.wavevectors3D","page":"Home","title":"CosmoMMF.wavevectors3D","text":"wavevectors3D([T=Float64], dims, box_size=(2π, 2π, 2π))\n\nArguments:\n\nT : output element type\ndims: array or tuple containing size\nbox_size: side lengths of box, normalizes the wavenumber\n\nReturns:\n\nTuple{Vector{T},Vector{T},Vector{T}}: kx, ky, kz wavevectors\n\nExamples\n\njulia> CosmoMMF.wavevectors3D((6,2,2))\n([0.0, 1.0, 2.0, -3.0, -2.0, -1.0], [0.0, -1.0], [0.0, -1.0])\n\n\n\n\n\n","category":"function"},{"location":"#CosmoMMF.xθ-Tuple{T} where T","page":"Home","title":"CosmoMMF.xθ","text":"xθ(x)\n\nHeaviside step times x. Returns x if the argument is positive, zero otherwise.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoMMF.θ-Tuple{T} where T","page":"Home","title":"CosmoMMF.θ","text":"θ(x)\n\nHeaviside step function. Returns one if the argument is positive, zero otherwise.\n\n\n\n\n\n","category":"method"}]
}
