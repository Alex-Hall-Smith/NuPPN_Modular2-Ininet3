##holds test functions for ppn and mppnp tests
import nugridse as mp
import ppn
import numpy as np

def dataget_ppn(case='', cycles=[]):
    '''
    Input
    ======
    case: string; name of ppn example being tested
    cycles: list; cycle numbers from ppn run to be tested (usually only a few)

    This function creates a test data directory and populates it with ppn data from the \
    CADC VOspace using the wget command. These data are the master results to be compared \
    against.

    '''
    import os

    os.mkdir('NuPPN/examples/regression_tests/' + case)

    for cycle in cycles:
        os.system("wget -q -P NuPPN/examples/regression_tests/" + case + " --content-disposition http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/projects/ppn/examples/" + case + "/iso_massf" + str(cycle) + ".DAT?view=data")


def dataget_mppnp(case='', files=[]):

    '''
     Input
    ======
    case: string; name of ppn example being tested
    files: list of file names to be downloaded
    data_type: either H5_out, H5_surf or H5_restart
    '''
    import os
    #
    path_case = 'NuPPN/examples/regression_tests/' + case

    os.mkdir(path_case)

    if 'out' in files[0]:
        mppnp_data_type = 'H5_out'
    elif 'surf' in files[0]:
        mppnp_data_type = 'H5_surf'
    if 'restart' in files[0]:
        mppnp_data_type = 'H5_restart'

    for file_name in files:
        os.system("wget -q -P " + path_case + " --content-disposition http://www.canfar.phys.uvic.ca/vospace/nodes/nugrid/data/projects/mppnp/examples/" + case + "/" + mppnp_data_type + "/"+file_name+"?view=data")

    print 'downloaded the following data into ',path_case,':'
    os.listdir(path_case)

def ppn_abundance_test(case='', cycles=[]):
    '''
    Docstring goes here
    '''

    p1=ppn.abu_vector(sldir='NuPPN/examples/regression_tests/' + case + '')
    p2=ppn.abu_vector(sldir='NuPPN/examples/' + case + '')

    for k in range(len(cycles)):
      cycle=cycles[k]
      b=p1.get('ABUNDANCE_MF', cycle)
      t=p2.get('ABUNDANCE_MF', cycle)

      if np.allclose(b, t, rtol=1e-5, atol=0, equal_nan=True):
         print "Test of iso_massf" + str(cycle) + ".DAT was successful."  
      else:
        print "Test of iso_massf" + str(cycle) + ".DAT failed."         

def mppnp_H5_test(cycles,case,files):

    '''
        Compares abundance profiles for cycles.
	...
    '''

    sefiles_ref =mp.se('NuPPN/examples/regression_tests/'+case)

    if 'out' in files[0]:
        mppnp_data_type = 'H5_out'
    elif 'surf' in files[0]:
        mppnp_data_type = 'H5_surf'
    if 'restart' in files[0]:
        mppnp_data_type = 'H5_restart'

    sefiles = mp.se('NuPPN/examples/'+case+'/'+mppnp_data_type+'/')


    for k in range(len(cycles)):
	cycle=cycles[k]
        t=sefiles.get(cycle,'iso_massf')
        b=sefiles_ref.get(cycle,'iso_massf')
        if np.allclose(b,t,rtol=0, atol= 1e-6, equal_nan=True):
	    if mppnp_data_type == 'H5_out' or mppnp_data_type == 'H5_restart':
            	print "Test of abundance profile of cycle ",cycle," was successful." 
	    else:
		print "Test of surface abundance of cycle ",cycle," was successful." 
        else:
	    if mppnp_data_type == 'H5_out' or mppnp_data_type == 'H5_restart':
            	print "Test of abundance profile of cycle ",cycle," FAILED." 
	    else:
		print "Test of surface abundance of cycle ",cycle," FAILED." 


