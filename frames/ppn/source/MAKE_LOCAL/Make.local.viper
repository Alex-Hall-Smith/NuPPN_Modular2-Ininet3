# Copy this file (or any of the machine specific Make.local.machines

# files) to Make.local and edit to suit your local configuration



# identifier for this local

LOCAL = viper

#Architecture 

ARCH =Linux_x86_64_gfortran


# make sure PHYSICS and SOLVER are absolute path names, for example

# specify the location of a full svn tree if you have one, like this:

PPN = /home/<UID>/NuPPN



# you don't need to change the following two but you can
# these are master settings -- modular2 doesn't have them
#PHYSICS = $(PPN)/physics/phys08
#SOLVER = $(PPN)/solver/NFR

PHYSICS = $(PPN)/physics
SOLVER = $(PPN)/solver


# where is the se library?

SEHOME = /home/<UID>/se



# where is mpi?

#MPIHOME = /here/is/your/openmpi



USE_SUPERLU = YES
# Location/version up to date as of 23/1/19
BLASLIB= -L/trinity/clustervision/CentOS/7/apps/OpenBLAS/0.2.18/lib/ -lopenblas



