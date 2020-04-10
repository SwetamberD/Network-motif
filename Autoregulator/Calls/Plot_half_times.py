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
# beta = 1.3
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
    Np_all = np.array([1, 2, 3])
    open_dir = "../Data/"
    for Np in Np_all:
        fname = "Half_times_data_for_Np_%d" % Np
        fname = fname.replace(".", "_")
        fname = fname + ".h5"
        hf = h5py.File(open_dir + fname, 'r')
        Nd_all = hf.get('Nd')
        Nd_all = np.array(Nd_all)
        Half_times = hf.get('Avg_half_times')
        Half_times = np.array(Half_times)
        hf.close()
        legend_text = r"$N_p = %d$" % Np
        plt.plot(Nd_all, Half_times, label=legend_text)
    plt.xlabel(r'$N_d$')
    plt.ylabel(r'Half times')
    plt.xlim([10**0, 5*10**4])
    plt.ylim([0.0, 80.0])
    plt.legend()
    plt.xscale('log')
    # plt.show()
    if save_fig:
        figname = "Half_times_vs_decoys_with_Nps_1_2_3.pdf"
        dirname = "../Plots/"
        plt.savefig(dirname + figname)
        print("saved:", figname)

    plt.show()

########################################################################################################################
def model_positive_feedback(vec, t, Np, Nd):

    x = vec[0]
    m = vec[1]

    # Parameters
    kappa1 = 0.25*100
    kappa2 = 1.3*100

    beta = 0.5

    r1 = 1.0
    r2 = 1.0#eps

    Np_r = Np*0.01
    Nd_r = Nd*0.01


    p_x  = 1.0 + 4.0*r1*x + 4.0*Np_r*x/(1+x**2)**2 + 4.0*Nd_r*r2*x/(1+r2*x**2)**2

    dx = -beta * (x - m)/p_x
    dm = kappa1 * Np_r/(1 + x** 2) + kappa2 * Np_r*x** 2/(1 + x** 2)  - m

    return [dx, dm]

########################################################################################################################
if __name__ == "__main__":
    main()