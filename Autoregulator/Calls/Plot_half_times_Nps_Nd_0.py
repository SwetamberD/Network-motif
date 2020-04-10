from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import h5py

########################################################################################################################
# Parameters as in MilJohWei2008 for positive feedback

# following parameters have the unit of min inverse
# alpha = 0.025
# beta = 1.7
# gamma_p = 0.05
# sigma = 0.5
# C = 10^-9 # mol inv
# gamma_m = 0.1

# following parameters are ratio (no units)
# cp = 10^7
# cd = 10^7


########################################################################################################################
def main():
    plot_half_times(save_fig = True)
    
########################################################################################################################
def plot_half_times(save_fig=True):
    open_dir = "../Data/"
    fname = "Half_times_for_different_Nps_Nd_0.h5"
    hf = h5py.File(open_dir+fname, 'r')
    Np_all = hf.get('Np')
    Np_all = np.array(Np_all)
    Half_times = hf.get('Half_times')
    Half_times = np.array(Half_times)
    hf.close()
    plt.plot(Np_all, Half_times,'o')
    plt.plot(Np_all, Half_times,'--k')
    plt.xlabel(r'$N_p$')
    plt.ylabel(r'Half times')
    if save_fig:
        figname = "Half_times_vs_decoys_with_Np_for_Nd_0.pdf"
        dirname = "../Plots/"
        plt.savefig(dirname + figname)
        print("saved:",figname)
    plt.show()
########################################################################################################################

########################################################################################################################
if __name__ == "__main__":
    main()