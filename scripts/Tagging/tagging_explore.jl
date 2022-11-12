
function get_clusbools(mass_of_average_cell,thresh, max_cats::AbstractArray{T,4}, den) where T
    """
    This function takes in a physical threshold float value and the 
    output of the maximum_signature function (4D-Array of signatures)
    to determine the cluster boolean filter.
    """
    
    # step 1. create a bool filter that is only true when 
    # cluster signatures are larger than threshold
    test_bools = (max_cats[:,:,:,1] .> thresh)
    
    nx, ny, nz = size(max_cats)[1:3]

    # step 2. tag components
    components = Images.label_components(test_bools)
    max_component = maximum(components)

    print("Components: $(max_component)\n")
    
    # step 3. loop through each component, assess average overdensity
    masses = zeros(T, max_component)
    volume = zeros(T, max_component)

    compare = zeros(nx,nx,nx);

    for i in 1:nx
        for j in 1:nx
            for k in 1:nx
                if components[i,j,k] != 0
                    masses[ components[i,j,k] ] += den[i,j,k]
                    volume[ components[i,j,k] ] += 1.0
                end
            end
        end
    end
    
    compbool = (masses .> 5e3 / mass_of_average_cell)  .& (volume .> 1)
    println( "viri frac: ", sum( (volume .> 1) .& (masses ./ volume .> 370) ) / max_component)
    println( "mass frac: ", sum(masses[compbool])   / sum(den) )
    
    clusbool = zeros(Bool, nx, nx, nx)
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                if components[i,j,k] != 0
                    if compbool[components[i,j,k]]
                        clusbool[i,j,k] = true
                    end
                end
            end
        end
    end
    
    clusbool
end

function load_clusbool(halo_kind, dm_only)
    if halo_kind == "Heavy"
        if dm_only == false
            @load "/global/cscratch1/sd/james12/halo_data/Heavy_Halos/"* "Heavy_Halos_Binned_Bool_TNG" * ".jld2" halo
        
        elseif dm_only == true
            @load "/global/cscratch1/sd/james12/halo_data/Heavy_Halos/"* "Heavy_Halos_Binned_Bool" * ".jld2" halo
        end
        
    elseif halo_kind == "Light"
        if dm_only == false
            @load "/global/cscratch1/sd/james12/halo_data/Light_Halos/"* "Light_Halos_Binned_Bool_TNG" * ".jld2" halo
        
        elseif dm_only == true
            @load "/global/cscratch1/sd/james12/halo_data/Light_Halos/"* "Light_Halos_Binned_Bool" * ".jld2" halo
        end
        
    elseif halo_kind == "UltraLight"
        if dm_only == false
            @load "/global/cscratch1/sd/james12/halo_data/Ultralight_Halos/final/"* "UltraLight_Halos_Binned_TNG_Bool" * ".jld2" halo
        
        elseif dm_only == true
            @load "/global/cscratch1/sd/james12/halo_data/Ultralight_Halos/final/"* "UltraLight_Halos_Binned_Bool" * ".jld2" halo
        end
        
    elseif halo_kind == "Both"
        if dm_only == false
            @load "/global/cscratch1/sd/james12/halo_data/"* "Both_Halos_Binned_Bool_TNG" * ".jld2" halo
        
        elseif dm_only == true
            @load "/global/cscratch1/sd/james12/halo_data/"* "Both_Halos_Binned_Bool" * ".jld2" halo
        end
        
    elseif halo_kind == "All"
        if dm_only == false
            @load "/global/cscratch1/sd/james12/halo_data/Heavy_p_Light_p_Ultralight/"* "All_Halos_Binned_TNG_Bool" * ".jld2" halo
        
        elseif dm_only == true
            @load "/global/cscratch1/sd/james12/halo_data/Heavy_p_Light_p_Ultralight/"* "All_Halos_Binned_Bool" * ".jld2" halo
        end
        
        
    end
    
    return halo
        
end

function get_sig_plot(sig_vec, f)
    """
    a function to calculate |dlog(M)/dlog(𝒮)| = ΔM^2 which is 
    known as the mass change with respect to signature. This can be 
    calculated for both filaments and clusters.
    """
    log10S = collect(-2:0.2:3) #this just makes an array [-2.0, -1.8, -1.6, ... , 2.8, 3.0]
    yy = zeros(Float64, size(log10S))
    for i in 1:size(yy,1)
        yy[i] = sum(
            f[sig_vec .> (10^log10S[i])]
            ) 
    end
    
    dydx = abs.( diff(yy.^2) ./ diff( log10S ) );
    midx = (log10S[1:size(log10S,1)-1] .+ log10S[2:size(log10S,1)]) ./ 2;
    return (10 .^ midx), dydx
end

function calc_structure_bools(combined_NEXUS::AbstractArray{T,4}, combined_NEXUSPLUS::AbstractArray{T,4}, den) where T
    """
    A function that returns the boolean filters for each structure based on physical criterion
    """
    
    N_cells = size(den,1) * size(den,2) * size(den,3)
    
    ####################################
    #     CLUSTER FILTER PROCESS
    ####################################
    
    halo_kind_val = "Light"
    dm_only_flag = true

    clusbool = load_clusbool(halo_kind_val, dm_only_flag);

    ####################################
    #    REMAINING FILTER PROCESS
    ####################################

    # we create 3 N^3 arrays with the signature corresponding to each structure
    cluster_signature = reshape(combined_NEXUS[:,:,:,1], N_cells)
    filament_signature = reshape(combined_NEXUSPLUS[:,:,:,2], N_cells)
    wall_signature = reshape(combined_NEXUSPLUS[:,:,:,3], N_cells);

    not_clus_flat = reshape( clusbool .== false, N_cells ); #bool array that is false where there are clusters
    filament_valid = filament_signature[ not_clus_flat ]; 
    flat_den_valid = reshape(den, N_cells)[ not_clus_flat ];

    S_fil, dM2_fil = get_sig_plot(filament_valid, flat_den_valid)
    max_dM, ind = findmax(dM2_fil) #find the peak of the mass change curve
    filament_thresh = S_fil[ind]

    wall_valid_filt = (not_clus_flat) .& (filament_signature .< filament_thresh)
    wall_valid = wall_signature[wall_valid_filt];
    wall_den_valid = reshape(den, N_cells)[wall_valid_filt];

    S_wall, dM2_wall = get_sig_plot(wall_valid, wall_den_valid)
    max_dM, ind = findmax(dM2_wall)
    wall_thresh = S_wall[ind]


    ####################################
    #        CREATE BOOL FILTERS
    ####################################

    cluster = den[clusbool]
    filbool = (combined_NEXUSPLUS[:,:,:,2] .> filament_thresh) .& (.!clusbool)
    filament = den[filbool]
    wallbool = (combined_NEXUSPLUS[:,:,:,3] .> wall_thresh) .& .!(filbool .| clusbool)
    wall = den[wallbool];
    
    
    return clusbool, filbool, wallbool, S_fil, dM2_fil, S_wall, dM2_wall
    
end
