import utils
import ppn
p=ppn.abu_vector('.')
import data_plot
from matplotlib.pylab import *

mp=p.get('mod')
sparse=10
cycles=mp[:1000:sparse]
form_str='%6.1F'
form_str1='%4.3F'

i=0
for cyc in cycles:
    T9  = p.get('t9',fname=cyc)
    Rho = p.get('rho',fname=cyc)
    mod = p.get('mod',fname=cyc)
#    time= p.get('agej',fname=cyc)*utils.constants.one_year
    time= p.get('agej',fname=cyc)
    close(i);figure(i);i += 1
    p.abu_chart(cyc,mass_range=[0,41],plotaxis=[-1,22,-1,22],lbound=(-6,0),show=False)
    title(str(mod)+' t='+form_str%time+'yr $T_9$='+form_str1%T9+' $\\rho$='+str(Rho))
    savefig('abu_chart_'+str(cyc).zfill(len(str(max(mp))))+'.png')
