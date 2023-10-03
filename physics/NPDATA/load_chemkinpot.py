from scipy.io import FortranFile


class chempot():
    def __init__(self):
        with open( 'chemkinpot.dat', 'r' ) as f:
            self.nt, self.nlrho, self.nye = np.fromfile(f, dtype=np.int32, count = 3)
            self.t9 = np.fromfile(f, dtype=np.float64, count = nt)
            self.ye = np.fromfile(f, dtype=np.float64, count = nye)
            self.lrho = np.fromfile(f, dtype=np.float64, count = nlrho)
            self.dt9, dlrho, dye = np.fromfile(f, dtype=np.float64, count = 3)
            self.mu_p = np.fromfile(f, dtype=np.float64, count = nt * nlrho * nye ).reshape([nt,nlrho,nye],order='F')
            self.mu_n = np.fromfile(f, dtype=np.float64, count = nt * nlrho * nye ).reshape([nt,nlrho,nye],order='F')
            self.iters = np.fromfile(f, dtype=np.int32, count = nt * nlrho * nye ).reshape([nt,nlrho,nye],order='F')
            self.ierr = np.fromfile(f, dtype=np.int32, count = nt * nlrho * nye ).reshape([nt,nlrho,nye],order='F')

        #with FortranFile( 'chemkinpot.dat', 'r' ) as f:
            #self.nt, self.lrho, self.nye = f.read_ints( dtype = np.int32 )
            #self.t9, self.ye, self.lrho = f.read_reals( dtype = np.float64 )
            #self.dt9, self.dlrho, self.dye = f.read_reals( dtype = np.float64 )
            #self.mu_p  = f.read_reals( dtype = np.float64 )
            #self.mu_n  = f.read_reals( dtype = np.float64 )
            #self.iters = f.read_ints( dtype = np.int32 )
            #self.ierr  = f.read_ints( dtype = np.int32 )

    def plot_mu(self):
        figure();
        pcolor(c.ye,c.t9,c.mu_n[:,self.nlrho//2,:])
        cbar = colorbar()
        cbar.set_label('$\mu_\mathrm{n}$')
        xlabel('$Y_\mathrm{e}$')
        ylabel('$T_9$')
        figure();
        pcolor(c.ye,c.t9,c.mu_p[:,self.nlrho//2,:])
        cbar = colorbar()
        cbar.set_label('$\mu_\mathrm{p}$')
        xlabel('$Y_\mathrm{e}$')
        ylabel('$T_9$')


c = chempot()
