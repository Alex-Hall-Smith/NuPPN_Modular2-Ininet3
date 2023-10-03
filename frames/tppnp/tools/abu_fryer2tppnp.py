#!/usr/bin/env python
import numpy as np
from sys import stdout

nzones = 1006
cfn = '/home/swj/work/data/tppnp/run_from_chris/forsydney/cfabundance.dat'
outfile = '/home/swj/work/data/tppnp/run_from_chris/forsydney/abundances.dat'

nams = ['NEUT ', 'H   1', 'HE  4', 'C  12', 'C  13', 'N  14', 'O  16',
        'NE 20', 'MG 24', 'SI 28', 'S  32', 'AR 36', 'CA 40', 'TI 44',
        'CR 48', 'FE 52', 'FE 54', 'NI 56', 'FE 56']

data = np.loadtxt(cfn)
nspecies = len(data[0])

with open(outfile, 'w') as f:
    for i in range(nzones):
        done = float(i) / float(nzones) * 100.
        done = int(done)
        stdout.write(" Converting abundance file. Progress: \
                      %s%%      %s"%(done,"\r"))
        stdout.flush();
        f.write('p' + str(i + 1) + '          ' + str(nspecies) + ' \n')
        for j in range(nspecies):
            f.write(nams[j] + "    {:.5e}".format(data[i][j]) + '\n')

