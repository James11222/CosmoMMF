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
# using CosmoMMF
include("tagging_explore.jl")

config = YAML.load_file("../../config.yaml")

############################################################
#                  Load Relevant Data
############################################################

print("Importing Density Field... \n")
#Density Field
@load config["input_directory"]*config["run_name"]*".jld2" den

print("Importing NEXUS Signatures... \n")
#Signature Arrays
if config["simulation_type"] == "Dark"
    
    save_name = config["run_name"] * "_NEXUS"
    @load config["output_directory"] * "DM_only_1024_snap_099/" * "max_sigs_" * save_name * "_full_signatures.jld2" sigmax
    sigmax_NEXUS = sigmax

    print("Importing NEXUS+ Signatures... \n")
    save_name = config["run_name"] * "_NEXUSPLUS"
    @load config["output_directory"] * "DM_only_1024_snap_099/" * "max_sigs_" * save_name * "_full_signatures.jld2" sigmax
    sigmax_NEXUSPLUS = sigmax
    
elseif config["simulation_type"] == "All"
    save_name = config["run_name"] * "_NEXUS"
    @load config["output_directory"] * "Full_Species_1024_snap_099/" * "max_sigs_" * save_name * "_full_signatures.jld2" sigmax
    sigmax_NEXUS = sigmax

    print("Importing NEXUS+ Signatures... \n")
    save_name = config["run_name"] * "_NEXUSPLUS"
    @load config["output_directory"] * "Full_Species_1024_snap_099/" * "max_sigs_" * save_name * "_full_signatures.jld2" sigmax
    sigmax_NEXUSPLUS = sigmax
end

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
clusbool, filbool, wallbool, S_fil, dM2_fil, S_wall, dM2_wall = calc_structure_bools(sigmax_NEXUS, sigmax_NEXUSPLUS, den)

####################################
#        SAVE BOOL FILTERS
####################################

print("Saving Structure Boolean Filters... \n")
@save config["output_directory"]*config["run_name"]*"_cluster_bool_filter.jld2" clusbool
@save config["output_directory"]*config["run_name"]*"_filament_bool_filter.jld2" filbool
@save config["output_directory"]*config["run_name"]*"_wall_bool_filter.jld2" wallbool


