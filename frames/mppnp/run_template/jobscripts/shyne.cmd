#BSUB -n 64                 # number of tasks in job         
#BSUB -R "span[ptile=64]&&span[hosts=1]"
#BSUB -R "select[hname!='shyne']&&select[hname!='compute011']"
#BSUB -J M20slu64           # job name
#BSUB -o $PWD/out           # output file name in which %J is replaced by the job ID
#BSUB -e $PWD/err           # error file name in which %J is replaced by the job ID

#run the executable
/opt/platform_mpi/bin/mpirun ./mppnp.exe
