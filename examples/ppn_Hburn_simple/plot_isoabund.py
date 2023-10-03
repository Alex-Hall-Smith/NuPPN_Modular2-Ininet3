import ppn
import utils
p=ppn.abu_vector('.')
import data_plot
from matplotlib.pylab import *

mp=p.get('mod')
sparse=10
cycles=mp[::sparse]
form_str='%6.1F'
form_str1='%4.3F'

i=0
for cyc in cycles:
    T9=p.get('t9',fname=cyc)
    Rho=p.get('rho',fname=cyc)
    mod=p.get('mod',fname=cyc)
#    time= p.get('agej',fname=cyc)*utils.constants.one_year
    time= p.get('agej',fname=cyc)
    close(i);figure(i);i += 1
    p.iso_abund(cyc,decayed=False,stable=False,show=False)
    title(str(mod)+' t='+form_str%time+'yr $T_9$='+form_str1%T9+' $\\rho$='+str(Rho))
    savefig('iso_abund_'+str(cyc).zfill(len(str(max(mp))))+'.png')
