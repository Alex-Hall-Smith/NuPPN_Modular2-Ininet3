import ppn
import utils
from matplotlib.pylab import *
symbs=utils.symbol_list('lines2')
x=ppn.xtime('.')
specs=['PROT','HE  4','C  12','N  14','O  16']
i=0
for spec in specs:
    x.plot('time',spec,logy=True,logx=True,shape=utils.linestyle(i)[0],show=False,title='')
    i += 1
ylim(-5,0.2)
legend(loc=0)
xlabel('$\log t / \mathrm{min}$')
ylabel('$\log X \mathrm{[mass fraction]}$')
savefig('abu_evolution.png')
