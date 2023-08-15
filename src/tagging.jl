"""
    make_the_clusbool() - Documentation
    
    A function that creates a Boolean Filter which selects only
    the clusters in a given cosmological 1+δ density cube. This function accomplishes
    this by finding a threshold value for 𝒮_cluster by looking at the
    change in virialization fraction as a function of 𝒮_cluster. 
    
    Note: Grid Resolution is important here. The resolution of data and max_sigs should
    be sufficiently high that a voxel's physical size is < 1 Mpc/h so clusters can be resolved.
    
    Arguments: 
    
    data - [3D Float Array] - data refers to the δ+1 data.
    max_sigs - [4D Float Array] - the maximum signatures array from CosmoMMF.maximum_signature()
    verbose - [Boolean] - a flag to allow the function to be more helpful and verbose.
    Δ - [Int] - The overdensity parameter threshold for determining virialization. This
    parameter comes from a paper by Gunn & Gott on spherical collapse. Commonly
    used values that are physically motivated can be 370, 200, or 500 (for R_200 or R_500).
    
    Returns:
    
    clusbool - [3D Boolean Array] - a boolean filter which can be used to isolate clusters
    in a cosmological volume.
    
"""
function make_the_clusbool(data, max_sigs, verbose, Δ) 
    
    # if verbose println("---------------------------------------") end
    # if verbose println("   Creating a Cluster Boolean Filter   ") end
    # if verbose println("---------------------------------------\n") end
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #         Step 1. Create the Initial Cluster Boolean Filter 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    # if verbose println("Step 1: Creating the Initial Cluster Boolean Filter...") end
    
    nx,ny,nz = size(data)
    total_volume = nx*ny*nz

    hist, bin_edges = np.histogram(max_sigs[:,:,:,1], bins=np.logspace(-5,2,nx))
    Smin = bin_edges[1:nx-1][argmax(hist)] #we make a cutoff at the peak of the distribution of signatures

    initial_clusbool = (max_sigs[:,:,:,1] .> Smin) 

    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #               Step 2. Identify Potential Candidates 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    # if verbose println("Step 2: Identifying Potential Candidates...") end
    
    # next we tag components with a component finder
    # a 3D array with groups of ID #'s for each component at each coordinate
    # in 1D: [0,0,0,0,0,1,1,1,0,0,0,0,2,2,0,0,0,0....]
    
    components = Images.label_components(initial_clusbool) 
    total_clusters_detected = maximum(components)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #                Step 3. Summarize Candidate Cluster Statistics 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # if verbose println("Step 3: Summarizing Candidate Cluster Statistics...") end
    
    #now we make 1D arrays with summary properties for each potential cluster
    
    density = zeros(total_clusters_detected)
    volume = zeros(total_clusters_detected)
    signatures = zeros(total_clusters_detected)

    for i in 1:nx
        for j in 1:nx
            for k in 1:nx
                if components[i,j,k] != 0 

                    signatures[ components[i,j,k] ] += max_sigs[i,j,k, 1]
                    density[ components[i,j,k] ] += data[i,j,k] #δ + 1 = ρ/ρ_bar
                    volume[ components[i,j,k] ] += 1.0

                end
            end
        end
    end
    
    #compute the average signature value of each cluster candidate
    ave_signatures = signatures ./ volume

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #     Step 4. Find the Signature Threshold via Virialization Curve 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    # if verbose println("Step 4: Finding the Signature Threshold via Virialization Curve...") end

    Ss = 10 .^ Array(range(round(log10(Smin)), 1, length=500))

    y = zeros(500)
    for (i,s_th) in enumerate(Ss)
        y[i] = sum((density .> Δ) .& (ave_signatures .> s_th)) / sum((ave_signatures .> s_th))
    end

    index, best_value = CosmoMMF.find_close_value_index(y, 0.5) 
    S_th = Ss[index]
    
    if verbose
        plt.figure()
        plt.plot(Ss, y, color="red", linewidth=3)
        plt.hlines(y=best_value, xmin=Smin, xmax = 1e1, color="Gray", linestyle="--", alpha=0.8)
        plt.vlines(x=S_th, ymin = 0, ymax = 1, color="Gray", linestyle="--", alpha=0.8)
        plt.xlabel("Average Cluster Signature")
        plt.ylabel("Fraction of Valid Clusters")
        plt.xlim(Smin, 1e1)
        plt.ylim(0,1)
        plt.semilogx()
        plt.show()
    end
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #                Step 5. Refine the Cluster Boolean Filter
    #                    to Meet Virialization Requirement
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    # if verbose println("Step 5: Refining the Cluster Boolean Filter to Meet Virialization Requirement...") end
    
    clusbool = max_sigs[:,:,:,1] .> S_th
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #               Step 6. Summary Statistics Report and Return 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    # if verbose println("Step 6: Summarizing Statistics Report...\n") end
    
    if verbose == true
        volume_fraction = sum(clusbool) / total_volume
        mass_fraction = sum(data[clusbool]) / sum(data)
        
        println("---------------------------------------")
        println("     Cluster Boolean Filter Report     ")
        println("---------------------------------------\n")
        
        println( "Signature Threshold of Clusters: $S_th")
        println( "Fraction of Valid Clusters: ", best_value) 
        println( "Volume Fraction of Clusters: ", volume_fraction)
        println( "Mass Fraction of Clusters: ", mass_fraction)
        println("")
    
    end

    return clusbool

end



"""
    calc_mass_change() - Documentation

    helper function to calculate |dM^2/dlog(𝒮)| = ΔM^2 which is 
    known as the mass change with respect to signature. This can be 
    calculated for both filaments and clusters.
    
    Smin, Smax are in log10 (default is -2, 3)
    
    Arguments:
    
    sig_vec - [Vector or 1D Float Array] - Valid Signatures (flattened)
    data_vec - [Vector or 1D Float Array] - Density Data (flattened)

    Returns:

    S, ΔM_2 - [Vector or 1D Float Array] - vectors for the mass change plot.
    
    
    
"""
function calc_mass_change(sig_vec, data_vec, Smin, Smax)
    
    #Initialize our arrays
    #this just makes an array [-2.0, -1.9, -1.8, ... , 2.9, 3.0]
    log10S = collect(Smin:0.1:Smax) 
    M = zeros(Float64, size(log10S))
    
    #this sums up all the mass in a structure type as a function of log(𝒮)
    for i in 1:size(M,1)
        filter = sig_vec .> (10^log10S[i])
        M[i] = sum(data_vec[filter]) 
    end
    
    #compute the derivative |dM^2/dlog(𝒮)|
    ΔM_2 = abs.( diff(M.^2) ./ diff( log10S ) );
    
    #compute log(𝒮) array used in derivatives (midpoints)
    midx = (log10S[1:size(log10S,1)-1] .+ log10S[2:size(log10S,1)]) ./ 2;
    S = (10 .^ midx)
    
    return S, ΔM_2
    
end

"""
    calc_structure_bools() - Documentation

    A function that returns the boolean filters for each structure based on physical criterion. 
    This is an implementation of the NEXUS+ method found in Cautun et al. 2012. This algorithm works 
    hierarchically meaning clusters are identified first, then filaments, then walls, and lastly voids.
    The user can opt to use their own "cluster" boolean filter, in Sunseri et al. 2023 we used a halo
    catalogue instead of a cluster catalog to study the effects of baryons beyond the dark matter halos. 

    If the user doesn't have a prior/external cluster boolean filter, the function will construct one
    from the signatures array produced by CosmoMMF.maximum_signature() using the NEXUS+ algorithm. This
    algorithm finds a signature threshold by using the response of virial fraction to cluster signature 
    strength. 

    Once the cluster boolean filter is created, the rest of the algorithm uses the mass change in filaments
    and walls as a function of signature to determine the threshold signature values for physical structures.
    If the user runs this function with verbose = true, they can also see these mass change curves.
    
    Arguments:
    
    data - [3D Float Array] - data refers to the δ+1 data.
    max_sigs - [4D Float Array] - the maximum signatures array
    verbose - [Boolean] - a flag for making the function be more verbose
    clusbool - [3D Boolean Array/nothing] - a premade boolean filter for clusters/halos. 
    This can be made externally or choose nothing if you want the implemented 
    NEXUS+ based cluster finder.
    Smin, Smax - Float - are in log10 (default is -3, 2)
    Δ - [Int] - The overdensity parameter threshold for determining virialization. This
    parameter comes from a paper by Gunn & Gott on spherical collapse. Commonly
    used values that are physically motivated can be 370, 200, or 500 (for R_200 or R_500).
    
    Returns:
    
    clusbool - [3D Boolean Array] - a boolean filter which isolates clusters in a cosmological volume.
    filbool - [3D Boolean Array] - a boolean filter which isolates filaments in a cosmological volume.
    wallbool - [3D Boolean Array] - a boolean filter which isolates walls in a cosmological volume.
    voidbool - [3D Boolean Array] - a boolean filter which isolates voids in a cosmological volume.

    if verbose = true
        S_fil, dM2_fil - [1D Float Arrays] - filament mass change curve variables for plotting
        S_wall, dM2_wall - [1D Float Arrays] - wall mass change curve variables for plotting
    
"""
function calc_structure_bools(data, max_sigs::AbstractArray{T,4}, verbose, clusbool=nothing, Smin=-3, Smax=2, Δ=370) where T
    
    N_cells = prod(size(data))
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #                 Step 1. Create Cluster Boolean Filter
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    if clusbool == nothing
        clusbool = make_the_clusbool(data, max_sigs, verbose, Δ)
    end

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #                 Step 2. Create Filament Boolean Filter
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # we create 3 N^3 arrays with the signature corresponding to each structure
    cluster_signature = reshape(max_sigs[:,:,:,1], N_cells)
    filament_signature = reshape(max_sigs[:,:,:,2], N_cells)
    wall_signature = reshape(max_sigs[:,:,:,3], N_cells);
    
    #Isolate Valid Filaments (not clusters)
    not_clus_flat = reshape( clusbool .== false, N_cells ); #bool array that is false where there are clusters
    filament_valid = filament_signature[ not_clus_flat ]; 
    flat_data_valid = reshape(data, N_cells)[ not_clus_flat ];

    #Compute Mass Change Curves and Find Filament Threshold
    S_fil, dM2_fil = calc_mass_change(filament_valid, flat_data_valid, Smin, Smax)
    max_dM, ind = findmax(dM2_fil) 
    filament_thresh = S_fil[ind]
    filbool = (max_sigs[:,:,:,2] .> filament_thresh) .& (.!clusbool)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #                   Step 3. Create Wall Boolean Filter
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    #Isolate Valid Walls (not clusters or filaments)
    wall_valid_filt = (not_clus_flat) .& (filament_signature .< filament_thresh)
    wall_valid = wall_signature[wall_valid_filt];
    wall_data_valid = reshape(data, N_cells)[wall_valid_filt];

    #Compute Mass Change Curves and Find Wall Threshold
    S_wall, dM2_wall = calc_mass_change(wall_valid, wall_data_valid, Smin, Smax)
    max_dM, ind = findmax(dM2_wall)
    wall_thresh = S_wall[ind]
    wallbool = (max_sigs[:,:,:,3] .> wall_thresh) .& .!(filbool .| clusbool)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #                   Step 4. Create Void Boolean Filter
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    voidbool = (clusbool .+ filbool .+ wallbool) .== 0
    
    
    if verbose
        
        println("---------------------------------------")
        println("     Filament Boolean Filter Report    ")
        println("---------------------------------------\n")
        
        volume_fraction = sum(filbool) / N_cells
        mass_fraction = sum(data[filbool]) / sum(data)
        
        println("Signature Threshold of Filaments: $filament_thresh")
        println( "Volume Fraction of Filaments: ", volume_fraction)
        println( "Mass Fraction of Filaments: ", mass_fraction)
        println("")


        println("---------------------------------------")
        println("       Wall Boolean Filter Report      ")
        println("---------------------------------------\n")
        
        volume_fraction = sum(wallbool) / N_cells
        mass_fraction = sum(data[wallbool]) / sum(data)
    
        println("Signature Threshold of Walls: $wall_thresh")
        println( "Volume Fraction of Walls: ", volume_fraction)
        println( "Mass Fraction of Walls: ", mass_fraction)  
        println("")
        
        println("---------------------------------------")
        println("       Void Boolean Filter Report      ")
        println("---------------------------------------\n")
        
        volume_fraction = sum(voidbool) / N_cells
        mass_fraction = sum(data[voidbool]) / sum(data)
    
        println( "Volume Fraction of Voids: ", volume_fraction)
        println( "Mass Fraction of Voids: ", mass_fraction) 
        println("\n---------------------------------------\n")
        
        
        # Plot the Mass Change Curves
        plt.figure()
        
        plt.plot(S_fil, dM2_fil / maximum(dM2_fil), color="blue", 
            linewidth=3, label="Filament")
        plt.plot(S_fil, dM2_wall / maximum(dM2_wall), color="green", 
            linewidth=3, label="Wall")
        
        plt.vlines(x=filament_thresh, ymin = 0, ymax = 1, 
            color="blue", linestyle="--", linewidth=3)
        plt.vlines(x=wall_thresh, ymin = 0, ymax = 1, 
            color="green", linestyle="--", linewidth=3)
        
        plt.xlabel("Signature Strength")
        plt.ylabel("ΔΜ^2 (arbitrary units)")
        plt.semilogx()
        plt.legend()
        plt.show()
        
        return clusbool, filbool, wallbool, voidbool, S_fil, dM2_fil, S_wall, dM2_wall
    
    else
        
        return clusbool, filbool, wallbool, voidbool
        
    end
    
end