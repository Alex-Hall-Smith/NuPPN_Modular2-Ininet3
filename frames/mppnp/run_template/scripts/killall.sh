# killall.sh script
# part of the mppnp utilities
# usage:
#   in the RUN directory use the run script as an
#   argument, e.g.
#   $ killall.sh run.helix_fa
# This will findout what the executable name is, and on 
# which machines you are running the job. The script will 
# then issue a killal command via ssh to all machines to 
# make sure the job properly terminates on all clients.
#
# Falk Herwig/NuGrid, April 2010

RUNSCRIPT=$1
NAME=`grep EXE= $RUNSCRIPT`
EXE=`echo ${NAME:4}`
echo The executable is $EXE ...

NAMEM=`grep MACHINEFILE= $RUNSCRIPT|grep -v NONE`
MACHINES=`echo ${NAMEM:12}`
echo The machines file is $MACHINES ...

n=0
while [ $n -lt $(wc -l <$MACHINES) ]
do
  let n=n+1
  set -- `head -n $n $MACHINES | tail -n 1`
  machine=$1
  echo Killing all jobs $EXE on machine $machine ....
  if [ `hostname` = $machine.phys.UVic.CA ] 
  then
      killall $EXE
  else
      ssh $machine killall $EXE
  fi
done
