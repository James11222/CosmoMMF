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

config = YAML.load_file("config.yaml")

############################################################
#                  Load Relevant Data
############################################################

#Density Field
@load config["output_directory"]*config["run_name"]*"_"*string(config["snapfile_root"])*".jld2" den

den = den[1:256,1:256,1:256]

saved_name = config["run_name"]*"_"*string(config["snapfile_root"])
#Signature Arrays
@load config["output_directory"]*"Sigs_Nexus_"*saved_name*".jld2" sig_array
combined_NEXUS = sig_array
@load config["output_directory"]*"Sigs_NexusPlus_"*saved_name*".jld2" sig_array
combined_NEXUSPLUS = sig_array

############################################################
# Simulation specifications (From Illustris website)
############################################################

# TNG300-3-Dark specifications
DM_particle_mass = 0.302538487429177 # in units of 1e10 Msun/h
N_DM = 244140625 
N_cells = size(den,1) * size(den,2) * size(den,3)
# an average grid cell has 
mass_of_average_cell = DM_particle_mass * N_DM / N_cells

# # TNG300-3 specifications
# DM_particle_mass = 0.302538487429177 # in units of 1e10 Msun/h
# N_DM = 244140625
# N_cells = size(den,1) * size(den,2) * size(den,3)
# # an average grid cell has 
# mass_of_average_cell = DM_particle_mass * N_DM / N_cells


####################################
#        CREATE BOOL FILTERS
####################################

clusbool, filbool, wallbool, S_fil, dM2_fil, S_wall, dM2_wall = CosmoMMF.calc_structure_bools(240.0, combined_NEXUS, combined_NEXUSPLUS,den)

####################################
#        SAVE BOOL FILTERS
####################################

@save config["output_directory"]*config["run_name"]*"_"*string(config["snapfile_root"])*"cluster_bool_filter.jld2" clusbool
@save config["output_directory"]*config["run_name"]*"_"*string(config["snapfile_root"])*"filament_bool_filter.jld2" filbool
@save config["output_directory"]*config["run_name"]*"_"*string(config["snapfile_root"])*"wall_bool_filter.jld2" wallbool

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

