#!/usr/bin/python
from __future__ import print_function
import os
import sdfpy as sp
import h5py
import numpy as np
from tqdm import tqdm

def snsph_sdf_to_hdf5(data_dir,prefix,files,output_filename):

	# read in list of zero padded files 
        with open(files) as f:
                content = f.readlines()
        files = [x.strip() for x in content]

	# get total number of files
        N_files = len(files)

	# make list of int from files
        for i in range(N_files):
        	files[i] = files[i].split('_')[1]

	files = np.asarray(files)

	# get information from last dump in list
	fin_data  = sp.SDFRead(str(data_dir)+str(prefix)+str(files[-1]))
	n_part = fin_data.parameters["npart"]
	tps = np.zeros((n_part,),dtype=np.int32)
	tps[:] = fin_data["ident"][:]
	pos_fin = np.zeros((n_part,3),dtype='float32')
        pos_fin[:,0],pos_fin[:,1],pos_fin[:,2] = fin_data["x"][:],fin_data["y"][:],fin_data["z"][:]
        tps_mass = np.zeros_like((tps),dtype='float32')
        tps_mass[:] = fin_data["mass"][:]
	
	# get final  positions
        ini_data = sp.SDFRead(str(data_dir)+str(prefix)+str(files[0]))
        idx = np.in1d(ini_data['ident'],tps)
	pos_ini = np.zeros((n_part,3),dtype='float32')
	pos_ini[:,0],pos_ini[:,1],pos_ini[:,2] = ini_data["x"][idx],ini_data["y"][idx],ini_data["z"][idx]

	# intialize our arrays
	rho = np.zeros((n_part,N_files),dtype='float32')
	temp = np.zeros((n_part,N_files),dtype='float32')
	time = np.zeros((N_files,),dtype='float32')

	# loop over all files in list only obtaining information
	# about the particles that survive to the end
        for j in tqdm(range(N_files)):
                
		# read in file
                sdf = sp.SDFRead(str(data_dir)+str(prefix)+str(files[j]))
        	idx = np.in1d(sdf['ident'],tps)

                # get global sim time to be written with each particle per file
                timeCF = sdf.parameters["timeCF"]
                elapsed_time = sdf.parameters["tpos"] * timeCF
		time[j] = elapsed_time
		rho[:,j] = sdf["rho"][idx]
		temp[:,j] = sdf["temp"][idx]

        # create hdf5 file
	hf = h5py.File(output_filename, 'w-')
	hf.create_dataset('time', data=time)
	hf.create_dataset('particle_mass', data=tps_mass)
        hf.create_dataset('pos_ini', data=pos_ini)
        hf.create_dataset('pos_fin', data=pos_fin)
	hf.create_dataset('particle_id', data=tps)
        hf.create_dataset('rho', data=rho)
	hf.create_dataset('temp', data=temp)
	hf.close()



# this program takes an arbitrary list of files of 
# snsph outout data and renames the list of files
# using padding for ease of use in reading in ordered 
def rename_snsph_dumps(data_dir,list_of_files,sim_name):

	with open(list_of_files) as f:
		content = f.readlines()
	files_orig = [x.strip() for x in content]
	files_new = files_orig[:]
	N_files = len(files_orig)
	
	for i in range(N_files):
		files_new[i] = files_new[i].split('.')[1]
		if len(files_new[i]) == 4 and float(files_new[i]) < 10000:
			files_new[i] = str(files_new[i]).zfill(5)
	
	files_new = np.asarray(files_new)
	
	for i in range(N_files):
		os.rename(str(data_dir)+str(files_orig[i]), str(data_dir)+str(sim_name)+'_'+str(files_new[i]))
	
	return


data_dir = '/lustre/scratch1/turquoise/carlnotsagan/Data/run15f1/'
output_filename = 'run15f1.h5'
files = './data_list.txt'
prefix = 'run15f1_'
sim_name = 'run15f1'
#rename_snsph_dumps(data_dir,files,sim_name)
snsph_sdf_to_hdf5(data_dir,prefix,files,output_filename)
