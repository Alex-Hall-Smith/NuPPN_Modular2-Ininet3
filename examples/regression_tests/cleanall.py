## Remove regression suite example directories
# Should we also clean run examples?

from os import listdir
from os.path import isdir
import shutil

path='NuPPN/examples/regression_tests/'

elements = listdir(path)

for element in elements:
    if isdir(path+element):
	#print 'delete ',path+element
	shutil.rmtree(path+element)         
