import illustris_python as il
import numpy as np
from tqdm import tqdm
import time 
import sys

#-------------------------------------------
#Initialize Data Directory
#-------------------------------------------
basePath = "/global/cscratch1/sd/james12/Illustris_TNG_Data/Dark_Sims/outputs" 

#-------------------------------------------
#Determine the Number of Halos 
#-------------------------------------------
halos = il.groupcat.loadHalos(basePath, snapNum=99, fields=None)
N_Halos = halos['count']

#-------------------------------------------
#Resolution & Constants of our 3D Histogram
#-------------------------------------------
L_Box = 205000 # check TNG300 - 1 for this
n = 1024 #make this small for debugging

pixel_size = 1 / 1024
m_dm_particle = 0.0047271638660809 * 1e10 # M_sun/h # for dm only


#-------------------------------------------
#Loop over each halo and add to the total
#-------------------------------------------

def bin_halos(weight_class):
    
    iters = 0
    
    if weight_class == "heavy":
        
        N_heavy = 0
        cube_heavy = np.zeros((1024,1024,1024))
        
        for i in tqdm(range(0, 5000)): #we know the first ~4000 halos are heavy
            dm_halo = il.snapshot.loadHalo(basePath, snapNum=99, id=i, partType="dm", fields='Coordinates')
            binned = np.floor((dm_halo / 205000) / pixel_size).astype(int)

            N_particles = binned.shape[0]
            M_halo = m_dm_particle * N_particles

            if M_halo >= 1e13:
                for j in range(0, N_particles):
                    x,y,z = binned[j]
                    cube_heavy[x,y,z] += 1
                N_heavy += 1

            
            iters += 1

            if N_heavy != iters:
                print(f"at halo {i} we are no longer loading in useful data...")
                break
                
        np.save("Heavy_Halos_Binned.npy", cube_heavy)
        
    elif weight_class == "light":
        
        N_light = 0
        cube_light = np.zeros((1024,1024,1024))
        
        for i in tqdm(range(4130, 40000)): #we know the first ~4000 halos are heavy
            dm_halo = il.snapshot.loadHalo(basePath, snapNum=99, id=i, partType="dm", fields='Coordinates')
            binned = np.floor((dm_halo / 205000) / pixel_size).astype(int)

            N_particles = binned.shape[0]
            M_halo = m_dm_particle * N_particles

            if M_halo < 1e13 and M_halo >= 1e12:
                for j in range(0, N_particles):
                    x,y,z = binned[j]
                    cube_light[x,y,z] += 1
                N_light += 1
            
            iters += 1

            if N_light != iters:
                print(f"at halo {i} we are no longer loading in useful data...")
                break
        
        np.save("Light_Halos_Binned.npy", cube_light)
        
    elif weight_class == "ultralight":
        
        N_ultralight = 0
        cube_ultralight = np.zeros((1024,1024,1024))
        
        for i in tqdm(range(35000, 400000)): 
            try:
                dm_halo = il.snapshot.loadHalo(basePath, snapNum=99, id=i, partType="dm", fields='Coordinates')
                binned = np.floor((dm_halo / 205000) / pixel_size).astype(int)

                N_particles = binned.shape[0]
                M_halo = m_dm_particle * N_particles
            except:
                N_particles = 0
                M_halo = 0

            if M_halo < 1e12 and M_halo > 1e11:
                for j in range(0, N_particles):
                    x,y,z = binned[j]
                    cube_ultralight[x,y,z] += 1
                N_ultralight += 1
                
            # if iters % 2000000 == 0:
            #     np.save(f"UltraLight_Halos_Binned_up_to_{iters}_iters.npy", cube_ultralight)
            
            iters += 1

        print(f"Number of Ultra Light Halos: {N_ultralight}")
        np.save("/global/cscratch1/sd/james12/Ultralight_Halos/final/UltraLight_Halos_Binned.npy", cube_ultralight)
        
    else:
        raise AssertionError("Invalid Weight Class")
                             

#-------------------------------------------
#Save Resulting histogram
#-------------------------------------------    

s = time.time()

weight_class_arg = sys.argv[1]
bin_halos(weight_class_arg)

sf = time.time()

print(f"Took {sf - s} Seconds")






