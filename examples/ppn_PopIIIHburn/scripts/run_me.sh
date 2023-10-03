# set the following parameters:
PPN_DIR=../..

function check_okay {
        if [ $error -ne 0 ]
        then
             echo $fail_warning
	     exit 1
        fi
}


[ -L $PPN_DIR/frames/ppn/NPDATA ]
error=$?; fail_warning="WARNING: $PPN_DIR/frames/ppn/NPDATA should be link and can not be found. Compile again in CODE."; check_okay

NPDATA_DIR=`readlink $PPN_DIR/frames/ppn/NPDATA`
echo "Using NPDATA in $NPDATA_DIR"

diff parameter.inc $PPN_DIR/frames/ppn/CODE/parameter.inc
error=$?
fail_warning="WARNING: parameter.inc is not the same as in CODE directory"
if [ $error -ne 0 ]
then
    echo $fail_warning
    echo "Attempting to fix this ..."
    cp parameter.inc $PPN_DIR/frames/ppn/CODE/
    cd $PPN_DIR/frames/ppn/CODE
    make distclean
    make
    cd -
    error=$?
    fail_warning="Error: Could not recompile automatically, try manually."
    check_okay
    echo "Automatic recompile apparently successful, now attempting to run ..."
fi

rm -f ../NPDATA
error=$?; fail_warning="WARNING: can not remove ../NPDATA"; check_okay
ln -s $NPDATA_DIR ..
error=$?; fail_warning="WARNING: can not link NPDATA"; check_okay

$PPN_DIR/frames/ppn/CODE/ppn.exe |tee fort.6

python plot_abu_evolution.py
python abu_chart.py 
python plot_isoabund.py

echo Finished making *.png plots. Compare with master_results/*png.

