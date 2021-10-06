########################################
#              Preamble
########################################
using PyCall
using FFTW
using LinearAlgebra
using JLD2
using YAML
using CosmoMMF

config = YAML.load_file("config.yaml")

########################################
#            Load in Data
########################################

@load config["output_directory"]*config["run_name"]*"_"*config["snapfile_root"]*".jld2" den

########################################
#  Run n' Save NEXUS/NEXUS+ Algorithm(s)
########################################

Rs = (√2) .^ (0:config["num_scales"])
print("\n loading different scales \n")

save_name = config["run_name"]*"_"*config["snapfile_root"]

if ARGS[1] == "NEXUS"
    @time combined_NEXUS = CosmoMMF.maximum_signature(Rs, den; alg=:NEXUS);
    CosmoMMF.save_max_sigs(config["output_directory"], "Sigs_Nexus_"*save_name, combined_NEXUS)
    combined_NEXUS = nothing
elseif ARGS[1] == "NEXUSPLUS"
    @time combined_NEXUSPLUS = CosmoMMF.maximum_signature(Rs, den; alg=:NEXUSPLUS);
    CosmoMMF.save_max_sigs(config["output_directory"], "Sigs_NexusPlus_"*save_name, combined_NEXUSPLUS)
    combined_NEXUSPLUS = nothing
elseif ARGS[1] == "BOTH"
    @time combined_NEXUS = CosmoMMF.maximum_signature(Rs, den; alg=:NEXUS);
    CosmoMMF.save_max_sigs(config["output_directory"], "Sigs_Nexus_"*save_name, combined_NEXUS)
    combined_NEXUS = nothing
    
    @time combined_NEXUSPLUS = CosmoMMF.maximum_signature(Rs, den; alg=:NEXUSPLUS);
    CosmoMMF.save_max_sigs(config["output_directory"], "Sigs_NexusPlus_"*save_name, combined_NEXUSPLUS)
    combined_NEXUSPLUS = nothing
else
    printstyled("Please give a valid algorithm input argument: \n NEXUS, NEXUSPLUS, or BOTH"; color = :red)
    
























########################################
#            Miscellaneous
########################################

# print(ARGS)

# x = parse(Int32,ARGS[3])
# y = parse(Int32,ARGS[4])

# print("\n")
# print(x+y)
# print("\n")