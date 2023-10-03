#!/bin/bash
#$ -S /bin/bash
# join stdout and stderr
#$ -j n
# change to working directory
#$ -cwd
# select parallel enviroment
#$ -pe mvapich48 192
# request resource, here the runtime
#$ -l h_rt=24:00:00
# bind the job to the requested cores
#$ -binding linear:48
#$ -q standard.q
#$ -m be
#$ -M samuel.jones@h-its.org
#$ -N tppnp
#$ -e .
#$ -o .
# load needed modules

cat $PE_HOSTFILE

source /etc/profile.d/modules.sh
module load tap
module load pso
module load sge
module load gcc/4.9.2
module load openmpi/gcc/64/1.6.3-qlc

mpirun -np ${NSLOTS} ./tppnp.exe
