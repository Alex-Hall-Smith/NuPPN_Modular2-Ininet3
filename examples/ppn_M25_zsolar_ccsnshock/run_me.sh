# select one of the 25 Msun cases from Raphael Hirschi, from the cases
# directory (just uncomment below):

case='he' # shocked He shell producing fe60
#case='c'  # shocked C shell producing fe60

#### most of the time no intervention required below this line #####

rm -f trajectory.input initial_abundance.dat
echo 'running case' $case 'using trajectory and initial abundance files:'
echo cases/$case/$case-trajectory.input

ln -s cases/$case/$case-trajectory.input  trajectory.input
ln -s cases/$case/$case-initial_abundance.dat  initial_abundance.dat

# run model with standard inlist
rm -f ppn_frame.input ppn_physics.input
ln -s ppn_frame_st.input ppn_frame.input
ln -s ppn_physics_st.input ppn_physics.input
./ppn.exe

# now decay
rm -f ppn_frame.input ppn_physics.input
ln -s cases/$case/$case-frame_dcy.input ppn_frame.input
ln -s ppn_physics_dcy.input ppn_physics.input
./ppn.exe
