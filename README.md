<p align="center">
  <img src="Images/CosmoMMF_Dark.png#gh-dark-mode-only" width="60%">
  <img src="Images/CosmoMMF_light.png#gh-light-mode-only" width="60%">
</p>



# CosmoMMF

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://James11222.github.io/CosmoMMF.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://James11222.github.io/CosmoMMF.jl/dev)
[![codecov](https://codecov.io/gh/James11222/CosmoMMF/branch/main/graph/badge.svg?token=cSdBAqnqya)](https://codecov.io/gh/James11222/CosmoMMF)
[![DOI](https://zenodo.org/badge/347520773.svg)](https://zenodo.org/badge/latestdoi/347520773)


Nexus Pipeline for analyzing the effects of baryonic matter on cosmological structures in IllustrisTNG simulations 

The `CosmoMMF.jl` package contains the algorithms necessary for a Multiscale Morphological Analysis (MMF) of cosmological simulations. The purpose of this package is to streamline our modified version of the NEXUS+ algorithm. We used this package in our work ([Sunseri et al. 2022](https://ui.adsabs.harvard.edu/abs/2023PhRvD.107b3514S/abstract)) to analyze the effects of baryonic matter on the Cosmic Web.

The NEXUS+ algorithm contains several steps as described in our paper ([Sunseri et al. 2022](https://ui.adsabs.harvard.edu/abs/2023PhRvD.107b3514S/abstract)). In general, we start with a density field, smooth it with a logarithmic Gaussian smoothing filter, then compute the hessian of the smoothed density field, use the eigenvalues of the hessian matrix to calculate the structure type signatures, find the maximum signatures over a range of smoothing scales, and apply physically based threshold criterion to categorize structures within the Cosmic Web. The entire package is implemented in `julia` and all of these steps are summarized inside of two functions. The first function `maximum_signature()` does the first several steps of the NEXUS+ algorithm to compute the maximum structure signatures, the second function is `calc_structure_bools()` which uses physical criteria to tag structures into 4 categories: clusters, filaments, walls, and voids. 

### General Code Usage

The general usage of the package would look like:

```julia
using CosmoMMF

@load "path/to/density_field.jld2" density_field; #load density field

Rs = (√2) .^ 0:10 #smoothing scales

max_signatures = CosmoMMF.maximum_signature(Rs, sphere_field, alg=:NEXUSPLUS) #compute maximum signatures

@load "path/to/cluster_boolean_filter.jld2" clusbool #load in externally computed boolean filter for clusters

clusbool, filbool, wallbool, S_fil, dM2_fil, S_wall, dM2_wall = CosmoMMF.calc_structure_bools(
                                                      clusbool, max_signatures, density_field) #tag structures
```

The output of `maximum_signature()` is a 4D Array where the 4th index denotes the signature type: 1 = clusters, 2 = filaments, 3 = walls. An example output of this can be seen below

<p align="center">
  <img src="Images/final_NEXUSPLUS_Signatures_hydro_dark.png#gh-dark-mode-only" width="100%">
  <img src="Images/final_NEXUSPLUS_Signatures_hydro.png#gh-light-mode-only" width="100%">
</p>

The boolean filters for each structure type produced by `calc_structure_bools()` can be used to tag structures within a density field, the results of this can be seen below

<p align="center">
  <img src="Images/final_tagging_figure_dark.png#gh-dark-mode-only" width="100%">
  <img src="Images/final_tagging_figure.png#gh-light-mode-only" width="100%">
</p>


### Additional Code Information 

* For exceptionally large arrays where containing the entire computation in memory is an issue, we include a `reduce_RAM_maximum_signature()` function which exploits the symmetries of the hessian matrix to reduce the memory footprint of the NEXUS+ algorithm. This function is used in the following

```julia
save_name = "save_name_string"
output_directory = "path/to/offload/temporary/calculations/"

max_signatures = CosmoMMF.reduce_RAM_maximum_signature(
                          Rs, output_directory, save_name, density_field, alg=:NEXUSPLUS)
```

* We also note in the `calc_structure_bools()` function, one does not have to load in their own boolean filter for clusters `clusbool`, one can set the first argument to be `nothing` and the code will use the criterion that a cluster must be virialized to be considered a valid cluster. We found this method to not be as accurate as other dedicated cluster finding algorithms, but we include it's implementation for completeness. For more information on this method, see [Cautun et al. 2013](https://academic.oup.com/mnras/article/429/2/1286/1038906).

