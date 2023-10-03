#!/bin/bash
# this is a runscript for running mppnp in MPI parallel mode on a number
# of machines; do not change anything in this file, instead before you 
# execute this script check the configuration parameters and comments 
# in config.sh
# FH, 141129

function check_okay {
        if [ $error -ne 0 ]
        then
             echo $fail_warning
	     exit 1
        fi
}

function check_info {
        if [ $error -ne 0 ]
        then
             echo $fail_warning
        fi
}

[ -f config.sh ] && source config.sh
error=$?
fail_warning="ERROR: no file config.sh. Run compile.sh script."
check_okay

# check if CODE is compiled with same parameter.inc as specified 
# in this RUN directory
CODE_DIR=$PPN_DIR/$FRAME/CODE
diff parameter.inc $CODE_DIR/parameter.inc
error=$?
fail_warning="WARNING: parameter.inc is not the same as in CODE directory. Need to recompile! Use compile.sh script."
check_okay

# check if PPN_DIR in Make.local is same as in this RUN directory
set -- `grep "PPN =" $PPN_DIR/$FRAME/CODE/Make.local` 
[[ `readlink -e $PPN_DIR` ==  *`echo $3`* ]]
error=$?
fail_warning="ERROR: PPN_DIR in this RUN dir and in CODE dir are not the same."
check_okay

# check if NPDATA is setup correctly
rm -f ../NPDATA
error=$?
fail_warning="INFO: can not remove link ../NPDATA - it may not exist."
check_info
PPN_DIR=`readlink -e $PPN_DIR`
ln -s $PPN_DIR/$PHYSICS/NPDATA ..
error=$?; fail_warning="WARNING: can not link NPDATA"; check_okay

NPDATA_DIR=`readlink ../NPDATA`
echo "Using NPDATA in $NPDATA_DIR"

# check if USEEPP directory is available and correctly setup
if [ -L ../USEEPP ]; then 
  USEEPP_DIR=`readlink ../USEEPP`
  echo "USEEPP is a link to $USEEPP_DIR"
elif [ -d ../USEEPP ] ; then
  echo "Using USEEPP in ../USEEPP"
else
  ln -s $PPN_DIR/$FRAME/USEEPP ../USEEPP
  error=$?; fail_warning="WARNING: can not create ../USEEPP link"; check_okay
fi

# create output directories if needed
[ -d H5_surf ]    || mkdir H5_surf
[ -d H5_out ]     || mkdir H5_out
[ -d H5_restart ] || mkdir H5_restart

# launching run
echo "Starting run with mpi "$MPIEXE" and mppnp executable "$EXE 
if [ $MACHINEFILE = 'NONE' ]
then
    echo "Not using machine file ...."
    $MPIEXE -np $NPROC   nice -n $NVAL ./$EXE > out  2> err.log &
else
    echo "Using machine file ...."
    echo "Launching: $MPIEXE -n $NPROC --bynode  --hostfile $MACHINEFILE  nice -n $NVAL ./$EXE > out  2> err.log &"
    $MPIEXE -n $NPROC --bynode  --hostfile $MACHINEFILE  nice -n $NVAL ./$EXE > out  2> err.log &
fi

