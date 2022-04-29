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

dm_file = "/global/cscratch1/sd/james12/Illustris_TNG_Data/TNG_Density_Cubes/simple_density_cubes/density_cube_snap_099_full_1024_dm.npy"
dm_density_cube = np.load(dm_file)

#load in Dark Matter Structure bool filters
clusbool = np.load("/global/cscratch1/sd/james12/NEXUS_analysis/quartiles/snap_099_full_1024_dm_q4_bool_nodes.npy")
filbool = np.load("/global/cscratch1/sd/james12/NEXUS_analysis/quartiles/snap_099_full_1024_dm_q3_bool_filaments.npy")
wallbool = np.load("/global/cscratch1/sd/james12/NEXUS_analysis/quartiles/snap_099_full_1024_dm_q2_bool_walls.npy")

full_density_cube = np.load(config["input_directory"] + config["run_name"] + ".npy") #load in baryon sim data


if (dm_density_cube == full_density_cube).all():
    raise AssertionError("dm only and full are the same...")
    
def compute_power_spectra(incremental_flag):
    """
    A function that computes power spectra by replacing baryonic
    matter using boolean filters from the quartiles file. 
    Saves various cross power files.
    
    Args:
        incremental_flag: (bool) a flag that if true incrementally replaces components instead of just only replacing 
        one component at a time.
        
    Returns:
        None
    """

    print("creating swapped clusters... \n")
    filtered_clus = dm_density_cube.copy()
    filtered_clus[clusbool] = full_density_cube[clusbool]

    if (filtered_clus == dm_density_cube).all():
        raise AssertionError("filtered_clus == dm_density_cube...")  

    print("creating swapped filaments... \n")
    filtered_fil = dm_density_cube.copy()
    
    if incremental_flag == True:
        filtered_fil[clusbool] = full_density_cube[clusbool]
        
    filtered_fil[filbool] = full_density_cube[filbool]

    if (filtered_fil == dm_density_cube).all():
        raise AssertionError("filtered_fil == dm_density_cube...")

    print("creating swapped walls... \n")
    filtered_wall = dm_density_cube.copy()
    
    if incremental_flag == True:
        filtered_wall[clusbool] = full_density_cube[clusbool]
        filtered_wall[filbool] = full_density_cube[filbool]
        
    filtered_wall[wallbool] = full_density_cube[wallbool]

    if (filtered_clus == dm_density_cube).all():
        raise AssertionError("filtered_wall == dm_density_cube...")
    else:
        print("don't worry they are all different")


    #-----------------------------------------
    #        Create NbodyKit Meshes
    #-----------------------------------------

    print("Creating Nbodykit Meshes...")
    dm_density_mesh = nb.ArrayMesh(dm_density_cube, BoxSize=config["xmax"])
    clus_NEXUS_mesh = nb.ArrayMesh(filtered_clus, BoxSize=config["xmax"])
    fil_NEXUS_mesh = nb.ArrayMesh(filtered_fil, BoxSize=config["xmax"])
    wall_NEXUS_mesh = nb.ArrayMesh(filtered_wall, BoxSize=config["xmax"])
    full_density_mesh = nb.ArrayMesh(full_density_cube, BoxSize=config["xmax"])

    #-----------------------------------------
    #       Compute Power Spectra
    #-----------------------------------------

    print("Computing Power Spectra...")
    print(".     |\n")
    result_all_dm = nb.FFTPower(first=dm_density_mesh, second=dm_density_mesh, mode='1d')
    print("..    |\n")
    result_NEXUS_clus = nb.FFTPower(first=dm_density_mesh, second=clus_NEXUS_mesh, mode='1d')
    print("...   |\n")
    result_NEXUS_fil = nb.FFTPower(first=dm_density_mesh, second=fil_NEXUS_mesh, mode='1d')
    print("....  |\n")
    result_NEXUS_wall = nb.FFTPower(first=dm_density_mesh,second=wall_NEXUS_mesh, mode='1d')
    print("..... |\n")
    result_all_full = nb.FFTPower(first=dm_density_mesh,second=full_density_mesh, mode='1d')
    print("......|\n")

    print("Saving Results...")
    if incremental_flag == True:
        result_all_dm.save("inc_cross_power" + "_pk_all_dm.json")
        result_NEXUS_clus.save("inc_cross_power" + "_pk_clus.json")
        result_NEXUS_fil.save("inc_cross_power" + "_pk_fil.json")
        result_NEXUS_wall.save("inc_cross_power" + "_pk_wall.json")
        result_all_full.save("inc_cross_power" + "_pk_all_full.json")
        
    else:
        result_all_dm.save("iso_cross_power" + "_pk_all_dm.json")
        result_NEXUS_clus.save("iso_cross_power" + "_pk_clus.json")
        result_NEXUS_fil.save("iso_cross_power" + "_pk_fil.json")
        result_NEXUS_wall.save("iso_cross_power" + "_pk_wall.json")
        result_all_full.save("iso_cross_power" + "_pk_all_full.json")
        

    finish = time.time()

    print("Power Spectra Computed in {0:0.4} seconds".format(finish - start))
    
if __name__ == "__main__":
    compute_power_spectra(incremental_flag=False)
    compute_power_spectra(incremental_flag=True)

