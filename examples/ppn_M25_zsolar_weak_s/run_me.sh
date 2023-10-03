# select one of the 25 Msun cases from Raphael Hirschi, from the cases
# directory (just uncomment below):

#case='central' # central trajectory
case='mixed'  # He core through C shell

#### most of the time no intervention required below this line #####

rm -f trajectory.input iniab.dat
echo 'running case' $case 'using trajectory and initial abundance files:'
echo cases/$case/$case-trajectory.input

ln -s cases/$case/$case-trajectory.input  trajectory.input

./ppn.exe
