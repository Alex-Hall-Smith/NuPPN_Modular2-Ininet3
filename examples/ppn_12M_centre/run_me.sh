# select one of the 12 Msun cases from Jones+ (2013), ApJ, from the cases
# directory (just uncomment below):

#case='all'
#case='h1_burn'
#case='he4_burn'
#case='c12_burn'
#case='ne20_burn'
case='o16_burn'
#case='si_burn'
#case='collapse'

#### most of the time no intervention required below this line #####

rm -f trajectory.input iniab.dat
echo 'running case' $case 'using trajectory and initial abundance files:'
echo cases/$case/$case-trajectory.input
echo cases/$case/$case-iniab.dat

ln -s cases/$case/$case-trajectory.input  trajectory.input
ln -s cases/$case/$case-iniab.dat iniab.dat

cat splash.ppn

./ppn.exe
