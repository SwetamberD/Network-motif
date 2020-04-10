from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import socket
import time
import h5py
###############################################################################################
# Main
###############################################################################################
def main():
    amp_ratio_computation(save_data=True)
    plot_amp_variation()

###############################################################################################
# Amplitude ratio computations
###############################################################################################

def amp_ratio_computation(save_data=True):
    vec = [10.0, 0.0, 10.0, 0.0, 0.0, 0.0]
    t = arange(0, 1000, 0.01)
    eps = 1.0
    # without decoys for reference
    Nd = 0.0
    Np = 1.0
    sol1 = odeint(model_repressilator, vec, t, args=(Np, Nd, eps))
    Amp_lower_bound_o = min(sol1[50000:80000, 0])
    Amp_upper_bound_o = max(sol1[50000:80000, 0])
    Amp_o = Amp_upper_bound_o - Amp_lower_bound_o
    Np_all = np.array([1.0, 2.0, 3.0])
    for Np in Np_all:
        Nd_all = np.linspace(1, 100000, 10000)
        Amp_ratio = []
        for Nd in Nd_all:
            sol2 = odeint(model_repressilator, vec, t, args=(Np, Nd, eps))
            Amp_lower_bound = min(sol2[50000:80000, 0])
            Amp_upper_bound = max(sol2[50000:80000, 0])
            Amp_1 = abs(Amp_upper_bound - Amp_lower_bound)
            Amp_rat = Amp_1 / Amp_o
            if Amp_1 > .01:
                #     plt.plot(sol1[30000:50000, 0])
                #     plt.plot(sol2[30000:50000, 0])
                #     plt.show()
                Amp_ratio.append(Amp_rat)
            else:
                Amp_ratio.append(0.0)
        if min(Amp_ratio) == 0.0:
            print("steady state detected!")
        if save_data:
            fname = "Amp_ratio_vs_decoy_Np_%d_eps_%1.1f_with_Np_ref_1" % (Np, eps)
            fname = fname.replace(".", "_")
            fname = fname + ".h5"
            hf = h5py.File(dirname + fname, 'w')
            hf.create_dataset('Nd', data=Nd_all)
            hf.create_dataset('Amp_ratio', data=Amp_ratio)
            hf.close()
            # print("saved:", fname)

###############################################################################################
# Plot
###############################################################################################

def plot_amp_variation():
    open_dir = "../Data/"
    Np_all = np.array([1.0, 2.0, 3.0])
    eps = 1.0
    for Np in Np_all:
        fname = "Amp_ratio_vs_decoy_Np_%d_eps_%1.1f_with_Np_ref_1" % (Np, eps)
        fname = fname.replace(".", "_")
        fname = fname + ".h5"
        label = r"$N_p = %1d$" % Np
        hf = h5py.File(open_dir + fname, 'r')
        decoy_nr = hf.get('Nd_all')
        decoy_nr = np.array(decoy_nr)
        Amp_ratio = hf.get('Amp_ratio')
        Amp_ratio = np.array(Amp_ratio)
        hf.close()
        plt.plot(decoy_nr,Amp_ratio, label=label)
    plt.xlabel(r"$N_d$")
    plt.ylabel("Amplitude ratio")
    plt.xscale('log')
    Amp_ratio_ref = decoy_nr * 0.0 + 1.0
    plt.plot(decoy_nr, Amp_ratio_ref, '--k')
    plt.legend()
    if save_fig:
        figname = "Amplitude_variation_with_Np_1_2_3_eps_%1.1f_Np_ref_1"% eps
        figname = figname.replace(".", "_")
        figname = figname + ".pdf"
        save_dirname1 = "../Plots/"
        # save_dirname2 = "../../../Fig_files/Fig_Repressilator/"
        plt.savefig(save_dirname1 + figname)
        # plt.savefig(save_dirname2 + figname)
        print("saved:", figname)
    #plt.show()
###############################################################################################
# Model_repressilator
###############################################################################################

def model_repressilator(vec,t, Np, Nd, eps):

    x1 = vec[0]
    m1 = vec[1]
    x2 = vec[2]
    m2 = vec[3]
    x3 = vec[4]
    m3 = vec[5]
    kappa = 800
    beta = 40
    r1 = 1.0
    r2 = 1.0
    Np_r = Np * 10 ** -2
    Nd_r = Nd * 10 ** -2

    p_x1 = 1.0 + 4.0 * r1 * x1 + 4.0 * Np_r * x1 / (1 + x1 ** 2) ** 2 + 4.0 * Nd_r * r2 * x1 / (1 + r2 * x1 ** 2) ** 2
    p_x2 = 1.0 + 4.0 * r1 * x2 + 4.0 * Np_r * x2 / (1 + x2 ** 2) ** 2 + 4.0 * Nd_r * r2 * x2 / (1 + r2 * x2 ** 2) ** 2
    p_x3 = 1.0 + 4.0 * r1 * x3 + 4.0 * Np_r * x3 / (1 + x3 ** 2) ** 2 + 4.0 * Nd_r * r2 * x3 / (1 + r2 * x3 ** 2) ** 2
    dx1 = -beta * (x1 - m1) / p_x1
    dm1 = kappa * Np_r / (1 + x3 ** 2) - m1

    dx2 = -beta * (x2 - m2) / p_x2
    dm2 = kappa * Np_r / (1 + x1 ** 2) - m2

    dx3 = -beta * (x3 - m3) / p_x3
    dm3 = kappa * Np_r / (1 + x2 ** 2) - m3

    return [dx1, dm1, dx2, dm2, dx3, dm3]

###############################################################################################
###############################################################################################

if __name__ == "__main__":
    main()