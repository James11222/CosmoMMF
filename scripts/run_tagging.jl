 ########################################
#              Preamble
########################################

using YAML
using Images
using PyCall
using FFTW
using LinearAlgebra
using JLD2
using YAML
using CosmoMMF

config = YAML.load_file("../config.yaml")

############################################################
#                  Load Relevant Data
############################################################

#Density Field
@load config["input_directory"]*config["run_name"]*".jld2" den

saved_name = config["run_name"]
#Signature Arrays
@load config["output_directory"]*"Sigs_Nexus_"*saved_name*".jld2" sig_array
combined_NEXUS = sig_array
@load config["output_directory"]*"Sigs_NexusPlus_"*saved_name*".jld2" sig_array
combined_NEXUSPLUS = sig_array

############################################################
# Simulation specifications (From Illustris website)
############################################################

# TNG-300 simulation specifications -> find average cell mass
if config["simulation_type"] == "Dark"
    DM_particle_mass = config["DM_particle_mass_Dark"] # in units of 1e10 Msun/h
    N_DM = config["N_DM"] 
    N_cells = size(den,1) * size(den,2) * size(den,3)
    # an average grid cell has 
    mass_of_average_cell = DM_particle_mass * N_DM / N_cells
elseif config["simulation_type"] == "All"
    DM_particle_mass = config["DM_particle_mass_TNG"] # in units of 1e10 Msun/h
    GAS_particle_mass = config["GAS_particle_mass_TNG"] # in units of 1e10 Msun/h
    N_DM = config["N_DM"] 
    N_GAS = config["N_GAS"] 
    N_cells = size(den,1) * size(den,2) * size(den,3)
    # an average grid cell has 
    mass_of_average_cell = (DM_particle_mass + GAS_particle_mass) * (N_DM + N_GAS) / N_cells
end


####################################
#        CREATE BOOL FILTERS
####################################

clusbool, filbool, wallbool, S_fil, dM2_fil, S_wall, dM2_wall = CosmoMMF.calc_structure_bools(mass_of_average_cell, 240.0, combined_NEXUS, combined_NEXUSPLUS, den)

####################################
#        SAVE BOOL FILTERS
####################################

@save config["output_directory"]*config["run_name"]*"_cluster_bool_filter.jld2" clusbool
@save config["output_directory"]*config["run_name"]*"_filament_bool_filter.jld2" filbool
@save config["output_directory"]*config["run_name"]*"_wall_bool_filter.jld2" wallbool

####################################
#        Make Plots
####################################

if ARGS[1] == "-p" 
    plt.ioff()

    f, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6))

    ax1.plot(S_fil, dM2_fil ./ maximum(dM2_fil), "ro")
    ax1.plot(S_fil, dM2_fil ./ maximum(dM2_fil), "r-")
    ax1.axvline(filament_thresh, color="white", linestyle="--")
    ax1.set_title("Filaments")
    ax1.set_ylabel("ΔM^2")
    ax1.set_xlabel("𝒮")
    ax1.set_xscale("log")

    ax2.plot(S_wall, dM2_wall ./ maximum(dM2_wall), "ro")
    ax2.plot(S_wall, dM2_wall ./ maximum(dM2_wall), "r-")
    ax2.xscale("log")
    ax2.set_title("Walls")
    ax2.set_ylabel("ΔM^2")
    ax2.set_xlabel("𝒮")
    ax2.axvline(wall_thresh, color="white", linestyle="--")
    f.tight_layout()
    f.save_fig
else
    print("plotting disabled, to turn on use -p flag")

