# --------------------------------------------------
#                    Preamble
# --------------------------------------------------

print("Reading in Imports... \n")

import scipy.ndimage as ndimg
import numpy as np
import yaml
import io

# Read YAML file
with open("../config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)
    
# --------------------------------------------------

density_file_dm = '/global/cscratch1/sd/james12/Illustris_TNG_Data/TNG_Density_Cubes/simple_density_cubes/density_cube_snap_099_full_1024_dm.npy'
density_file_all_species = "/global/cscratch1/sd/james12/Illustris_TNG_Data/TNG_Density_Cubes/simple_density_cubes/density_cube_snap_099_full_1024_all_species.npy" 

save_dir = "/global/cscratch1/sd/james12/NEXUS_analysis/quartiles/"

# --------------------------------------------------

def calc_structure_bools(density_file, type_flag):
    
    """
    Function that takes in a density cube file and
    returns 4 boolean filters, one for each quartile.
    
    Returns:
        q1_bool: (ndarray, dtype=bool) for voids
        q2_bool: (ndarray, dtype=bool) for walls
        q3_bool: (ndarray, dtype=bool) for filaments
        q4_bool: (ndarray, dtype=bool) for nodes
    """
    
    # total_volume = config["xmax"] ** 3
    # volume_per_cell = total_volume / (config["nx"]**3)
    
    print("Loading in Density File... \n")
    
    density_cube = np.load(density_file)
    smoothed_density = ndimg.gaussian_filter(density_cube, sigma=1)
    # mass_cube = density_cube * volume_per_cell
    # smoothed_mass = ndimg.gaussian_filter(mass_cube, sigma=0.5)
    
    print("Computing Mass Quartiles... \n")

    q1, q2, q3 = np.quantile(smoothed_density, [0.77, 0.95, 0.998]) #voids, walls, filaments, nodes
    
    print("Creating Boolean Filters... \n")

    structure_class = np.zeros(shape=smoothed_density.shape)
    structure_class[np.where (smoothed_density<q1)] = 0
    structure_class[np.where ((smoothed_density<q2) & (smoothed_density>q1))] = 1
    structure_class[np.where ((smoothed_density<q3) & (smoothed_density>q2))] = 2
    structure_class[np.where (smoothed_density>q3)] = 3

    q1_bool = structure_class == 0
    q2_bool = structure_class == 1
    q3_bool = structure_class == 2
    q4_bool = structure_class == 3
    
    print("Saving Boolean Filters... \n")
    
    
    if type_flag == "dm":
        np.save(save_dir + "snap_099_full_1024_dm_q1_bool_voids.npy", q1_bool)
        np.save(save_dir + "snap_099_full_1024_dm_q2_bool_walls.npy", q2_bool)
        np.save(save_dir + "snap_099_full_1024_dm_q3_bool_filaments.npy", q3_bool)
        np.save(save_dir + "snap_099_full_1024_dm_q4_bool_nodes.npy", q4_bool)
        
    elif type_flag == "all":
        # np.save(save_dir + "snap_099_full_1024_all_species_q1_bool_nodes.npy", q1_bool)
        # np.save(save_dir + "snap_099_full_1024_all_species_q2_bool_voids.npy", q2_bool)
        # np.save(save_dir + "snap_099_full_1024_all_species_q3_bool_walls.npy", q3_bool)
        # np.save(save_dir + "snap_099_full_1024_all_species_q4_bool_filaments.npy", q4_bool)
        print('skip')
        
    
if __name__ == "__main__":
    calc_structure_bools(density_file_dm, "dm")
    # calc_structure_bools(density_file_all_species, "all")
    


