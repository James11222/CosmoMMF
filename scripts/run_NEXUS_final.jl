########################################
#              Preamble
########################################
using PyCall
using FFTW
using LinearAlgebra
using JLD2
using YAML
using CosmoMMF

config = YAML.load_file("../config.yaml")

########################################
#            Load in Data
########################################

@load config["input_directory"]*config["run_name"]*".jld2" den
print("\n Loading in Density Cube \n")

########################################
#  Run n' Save NEXUS/NEXUS+ Algorithm(s)
########################################

Rs = (√2) .^ (0:config["num_scales"])
print("\n loading different scales \n")
output_dir = config["output_directory"]*ARGS[1]*"/"

if ARGS[1] == "NEXUS"
    save_name = config["run_name"] * "_NEXUS"
    @time CosmoMMF.reduce_RAM_maximum_signature(Rs, output_dir, save_name, den; alg=:NEXUS);
elseif ARGS[1] == "NEXUSPLUS"
    save_name = config["run_name"] * "_NEXUSPLUS"
    @time CosmoMMF.reduce_RAM_maximum_signature(Rs, output_dir, save_name, den; alg=:NEXUSPLUS);
else
    printstyled("Please give a valid algorithm input argument: \n NEXUS or NEXUSPLUS"; color = :red)
end
    

