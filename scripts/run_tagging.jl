########################################
#              Preamble
########################################

print("Importing Packages... \n")
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

print("Importing Density Field... \n")
#Density Field
@load config["input_directory"]*config["run_name"]*".jld2" den

print("Importing NEXUS Signatures... \n")
#Signature Arrays
save_name = config["run_name"] * "_NEXUS"
@load config["output_directory"] * "max_sigs_" * save_name * "_full_signatures.jld2" sigmax
sigmax_NEXUS = sigmax

print("Importing NEXUS+ Signatures... \n")
save_name = config["run_name"] * "_NEXUSPLUS"
@load config["output_directory"] * "max_sigs_" * save_name * "_full_signatures.jld2" sigmax
sigmax_NEXUSPLUS = sigmax

#reduce RAM usage by removing sigmax from
sigmax = nothing

############################################################
# Simulation specifications (From Illustris website)
############################################################

print("Importing Simulation Specifications... \n")

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

print("Calculating Structure Boolean Filters... \n")
clusbool, filbool, wallbool, S_fil, dM2_fil, S_wall, dM2_wall = CosmoMMF.calc_structure_bools(mass_of_average_cell, 240.0, sigmax_NEXUS, sigmax_NEXUSPLUS, den)

####################################
#        SAVE BOOL FILTERS
####################################

print("Saving Structure Boolean Filters... \n")
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

