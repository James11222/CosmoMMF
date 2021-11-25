import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
import yaml
import pprint
from nbodykit.lab import *
from nbodykit import setup_logging, style

snap_num = 99
snap_str = str(99)

config = {'nx' : 1024, 
          'xmax' : 205000.0, 
          'xmin' : 0.0,
          'input_dir' : '/global/cscratch1/sd/james12/Illustris_TNG_Data/Full_Sims/snapdir_099/',
          'output_dir' : '/global/cscratch1/sd/james12/Illustris_TNG_Data/TNG_Density_Cubes/'}

# config = yaml.load(open("config.yaml"))
# pp = pprint.PrettyPrinter()
# pp.pprint(config)
boxsize = config["xmax"] - config['xmin']



input_dir_Both = '/tigress/zequnl/TNG/300-1/'

f_gas = HDFCatalog(f'{input_dir_Both}/*.hdf5', root='PartType0')
f_dm = HDFCatalog(f'{input_dir_Both}/*.hdf5', root='PartType1')
f_star = HDFCatalog(f'{input_dir_Both}/*.hdf5', root='PartType4')
f_bh = HDFCatalog(f'{input_dir_Both}/*.hdf5', root='PartType5')

# f_dm['Masses'] = 0.254939359915231
f_dm['Masses'] = 0.00398342749867548 #add mass field for dark matter particles 10^10 M_sun/h
# f_gas['Coordinates'] *= 1e-3
# f_dm['Coordinates'] *= 1e-3
# f_star['Coordinates'] *= 1e-3
# f_bh['Coordinates'] *= 1e-3

print('done coordinate transformation!')

combined = MultipleSpeciesCatalog(['gas', 'dm', 'star', 'bh'], f_gas, f_dm, f_star, f_bh)
# # the combined mesh, weighted by mass
mesh = combined.to_mesh(
    Nmesh=2048, BoxSize=boxsize*1e-3, compensated=True, resampler='tsc', interlaced=True, weight='Masses', position='Coordinates')

# print('computing the FFT!')

# # compute the power, specifying desired linear k-binning
# r_all = FFTPower(mesh, mode='1d')
# Pk_all = r_all.power

# r_all.save('/tigress/zequnl/all.json')

