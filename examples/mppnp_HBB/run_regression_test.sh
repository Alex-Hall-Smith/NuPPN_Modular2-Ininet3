### Part of regression test suite. Do not execute manually!

#set a  minimum of 20 cycles
restartnum=`awk '{print $1}' last_restart.out`
modstop=$((restartnum+20))
echo 'modstop set to '$modstop
sed -i '/modstop/s/.*/        modstop='$modstop'/' ppn_frame.input
echo 'start running example '`pwd`
mpirun --allow-run-as-root -np 2 ./mppnp.exe  > out  2> err.log &
pid=$!

# If this script is killed, kill mppnp
trap "kill $pid 2> /dev/null" EXIT

while kill -0 $pid 2> /dev/null; do
    # Do stuff
    echo 'MPPNP run in progress...\r'
    sleep 10
done
#print output after finishing
cat out
cat err.log

# Disable the trap on a normal exit.
trap - EXIT
