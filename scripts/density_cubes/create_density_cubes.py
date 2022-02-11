import sys
import numpy as np
from nbodykit.lab import *
import h5py

snap_num = int(sys.argv[1])
snapdir_str = 'snapdir_' + str(snap_num).zfill(3) + '/'

if sys.argv[2] == 'dark':
    config = {'nx' : 1024, 
              'xmax' : 205000.0, 
              'xmin' : 0.0,
              'input_dir' : '/global/cscratch1/sd/jialiu/NEXUS_TNG/L205n2500TNG_DM/' + snapdir_str,
              'output_dir' : '/global/cscratch1/sd/james12/Illustris_TNG_Data/TNG_Density_Cubes/'}

elif sys.argv[2] == 'all':  
    config = {'nx' : 1024, 
              'xmax' : 205000.0, 
              'xmin' : 0.0,
              'input_dir' : '/global/cscratch1/sd/jialiu/NEXUS_TNG/L205n2500TNG/' + snapdir_str,
              'output_dir' : '/global/cscratch1/sd/james12/Illustris_TNG_Data/TNG_Density_Cubes/'}
    
else:
    print("Invalid Input")

#note: if we need bigfiles use BigFileCatalog instead of HDFCatalog and
#the file extension is .big

def compute_snap(snapnum, run_type):
    
    snapnum_str = str(snapnum).zfill(3)
    input_files_string = config['input_dir'] + 'snap_' + snapnum_str + '.*.hdf5'
    output_file_string = config['output_dir'] + 'density_cube_snap_' + snapnum_str + '_full'
    
    if run_type=='dark':
        
        # DM
        f = HDFCatalog(input_files_string)
        f['Position'] = f['PartType1/Coordinates']
        m = f.to_mesh(config['nx'], config["xmax"] - config['xmin'])
        field = m.compute()
        np.save(output_file_string + '_' + str(config['nx']) + "_dm.npy", field)
        
        
    elif run_type=='all':
        
        f_gas = HDFCatalog(input_files_string, root='PartType0')
        f_dm = HDFCatalog(input_files_string, root='PartType1')
        f_star = HDFCatalog(input_files_string, root='PartType4')
        f_bh = HDFCatalog(input_files_string, root='PartType5')
        
        
        f_dm['Masses'] = 0.00398342749867548 #add mass field for dark matter particles 10^10 M_sun/h
        
        combined = MultipleSpeciesCatalog(['gas', 'dm', 'star', 'bh'], f_gas, f_dm, f_star, f_bh)
        # # the combined mesh, weighted by mass
        mesh = combined.to_mesh(Nmesh=config['nx'], BoxSize=config["xmax"] - config['xmin'], 
                                compensated=True, resampler='tsc', interlaced=True, 
                                weight='Masses', position='Coordinates')
        
        full_field = mesh.compute()
        np.save(output_file_string + '_' +  str(config['nx']) +  "_all_species.npy", full_field)
        

    else:
        raise ValueError("Please give a valid run_type ('dark' or 'all')")
        


print("""\nComputing your density cubes for snap """ + str(snap_num) + '\n')

compute_snap(snap_num, sys.argv[2])


