
#set a  minimum of 20 cycles
sed -i 's/modstop  = 79991/modstop  =  78011/' ./ppn_frame.input
cat ppn_frame.input
echo 'start example in '`pwd`
mpirun --allow-run-as-root -np 2 ./mppnp.exe  > out  2> err.log &
pid=$!
#watch 'cat out; cat err.log; cat summaryinfo.dat'

# If this script is killed, kill the `cp'.
trap "kill $pid 2> /dev/null" EXIT

# While copy is running...
while kill -0 $pid 2> /dev/null; do
    # Do stuff
    echo -ne 'MPPNP run in progress...\r'
    sleep 5
done
#output if successufull finished
cat out
cat err.log

# Disable the trap on a normal exit.
trap - EXIT
