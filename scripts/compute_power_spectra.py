##########################################
#               PREAMBLE
##########################################

import numpy as np 
import yaml
import io
import nbodykit.lab as nb
import time

start = time.time()
#-----------------------------------------
#            Read YAML file
#-----------------------------------------

with open("../config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)
    
print("Configuration File Settings are: \n", config)
    
#-----------------------------------------
#   Read in Filters + Density Field
#-----------------------------------------

print("Reading in Data...")
density_cube = np.load(config["input_directory"] + config["run_name"] + ".npy")

clusbool = np.load(config["output_directory"] + config["run_name"] + "_cluster_bool_filter.npy")
filbool = np.load(config["output_directory"] + config["run_name"] + "_filament_bool_filter.npy")
wallbool = np.load(config["output_directory"] + config["run_name"] + "_wall_bool_filter.npy")

filtered_clus = density_cube
filtered_clus = filtered_clus * clusbool.astype(int) #compute Pk of clusters only

filtered_fil = density_cube
filtered_fil = filtered_fil * filbool.astype(int)
filtered_fil = filtered_clus + filtered_fil #compute Pk of clusters + filaments

filtered_wall = density_cube
filtered_wall = filtered_wall * wallbool.astype(int)
filtered_wall = filtered_clus + filtered_fil + filtered_wall #compute Pk of clusters + filaments + walls



#-----------------------------------------
#        Create NbodyKit Meshes
#-----------------------------------------

print("Creating Nbodykit Meshes...")
density_mesh = nb.ArrayMesh(density_cube, BoxSize=config["xmax"])

clus_NEXUS_mesh = nb.ArrayMesh(filtered_clus, BoxSize=config["xmax"])
fil_NEXUS_mesh = nb.ArrayMesh(filtered_fil, BoxSize=config["xmax"])
wall_NEXUS_mesh = nb.ArrayMesh(filtered_wall, BoxSize=config["xmax"])

#-----------------------------------------
#       Compute Power Spectra
#-----------------------------------------

print("Computing Power Spectra...")
result_all = nb.FFTPower(density_mesh, mode='1d')
result_NEXUS_clus = nb.FFTPower(clus_NEXUS_mesh, mode='1d')
result_NEXUS_fil = nb.FFTPower(fil_NEXUS_mesh, mode='1d')
result_NEXUS_wall = nb.FFTPower(wall_NEXUS_mesh, mode='1d')

# pk_all = result_all.power
# pk_clus = result_NEXUS_clus.power
# pk_fil = result_NEXUS_fil.power
# pk_wall = result_NEXUS_wall.power


result_all.save(config["run_name"] + "_pk_all.json")
result_NEXUS_clus.save(config["run_name"] + "_pk_clus.json")
result_NEXUS_fil.save(config["run_name"] + "_pk_fil.json")
result_NEXUS_wall.save(config["run_name"] + "_pk_wall.json")
# print("Saving Power Spectra...")
# (pd.DataFrame({"k": pk_all["k"], "power":pk_all['power'], "shotnoise":pk_all.attrs['shotnoise']}).to_csv(config["run_name"] + "_pk_all.csv")
# print(pk_all.attrs['shotnoise'])
# (pd.DataFrame({"k": pk_clus["k"], "power":pk_clus['power'], "shotnoise":pk_clus.attrs['shotnoise']}).to_csv(config["run_name"] + "_pk_clus.csv")
# print(pk_clus.attrs['shotnoise'])
# (pd.DataFrame({"k": pk_fil["k"], "power":pk_fil['power'], "shotnoise":pk_fil.attrs['shotnoise']}).to_csv(config["run_name"] + "_pk_fil.csv")
# print(pk_fil.attrs['shotnoise'])
# (pd.DataFrame({"k": pk_wall["k"], "power":pk_wall['power'], "shotnoise":pk_wall.attrs['shotnoise']}).to_csv(config["run_name"] + "_pk_wall.csv")
# print(pk_wall.attrs['shotnoise'])

#-----------------------------------------
#          Plot Results
#-----------------------------------------

# print("Plotting...")
# plt.figure(figsize=(12,8))
# # print the shot noise subtracted P(k)
# plt.loglog((pk_all['k']) * 1e3, (pk_all['power'].real - pk_all.attrs['shotnoise'])*1e-9, label='All Data', color="black") #Illustris units are kpc we need to convert to Mpc
# plt.loglog((pk_clus['k']) * 1e3, (pk_clus['power'].real - pk_clus.attrs['shotnoise'])*1e-9, label='NEXUS+ Cluster', color="red", linestyle="--")
# plt.loglog((pk_fil['k']) * 1e3, (pk_fil['power'].real - pk_fil.attrs['shotnoise'])*1e-9, label='NEXUS+ Filament', color="blue", linestyle="-.")
# plt.loglog((pk_wall['k']) * 1e3, (pk_wall['power'].real - pk_wall.attrs['shotnoise'])*1e-9, label='NEXUS+ Wall', color="green", linestyle=":")


# plt.title("Power Spectra: DM Only")
# # format the axes
# plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
# plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
# plt.legend(fontsize=14)
# plt.savefig("NEXUS+_pk_full_species.png", format="png", dpi=300, bbox_inches="tight")
# plt.show()

finish = time.time()

print("Power Spectra Computed in {0:0.4} seconds".format(finish - start))