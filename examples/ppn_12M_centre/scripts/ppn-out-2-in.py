import ppn

'''

This script takes an abundance vector from the .. directory, whose number is
given by <fname>, and writes a ppn input file with the name <outfile>, which
can be used with iabuini=12 in ppn_frame.input.
If you change the physics (e.g. nuclear reaction rates), run a burning stage,
and then want to propagate through to the next burning stage, you should use
this script to write a new initial abundance for the next stage and copy it
into the cases folder in the appropriate place.

'''

fname=518
outfile='collapse-iniab.dat'

p=ppn.abu_vector('..')
abunds=p.get('ABUNDANCE_MF',fname=fname)
isos=p.isotopes
f=open(outfile,'w')
for i in range(len(isos)):
    iso = isos[i]
    nam = iso.split('-')
    if str(nam[1][-2:]).upper() != 'M1' and  \
       str(nam[1][0]) != '*' and \
       str(nam[1][0]) != 'g' :
        s = 'D ' + str(nam[0]).upper().ljust(3) + \
            '    ' + str(nam[1]).upper().rjust(3) + \
            '    '+"{:e}".format(abunds[i])
        f.write(s+" \n")
f.close()
