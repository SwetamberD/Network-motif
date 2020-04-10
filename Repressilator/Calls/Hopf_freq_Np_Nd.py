from __future__ import division, print_function
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import time
from sympy import *
from sympy import Matrix

##############################################################################################
# main
###############################################################################################
def main():
    plot_freq_ratio_with_Nd(eps=1.0, save_fig=True)

###############################################################################################
# Plot Frequency ratio vs number of decoys
###############################################################################################

def plot_freq_ratio_with_Nd(eps=10.1,save_fig=True):
    Gamma = 40.0
    Nd_all = np.linspace(0, 10**5, 10000)
    Np_all = np.array([1.0, 2.0, 3.0])
    kappa = 800
    for Np in Np_all:
        omega_ratio_all = []
        for Nd in Nd_all:
            omega = freq_ratio(Np = Np, Nd=Nd, beta=Gamma,eps=eps, kappa=kappa)
            omega_ratio_all.append(omega)
        label = r"$N_p =%d$" % Np
        plt.plot(Nd_all, omega_ratio_all, label=label)
        plt.legend(loc="upper right")
    plt.xlabel(r"$N_d$")
    plt.ylabel(r"$\omega_r$")
    omega_ratio_ref = np.linspace(1,10**6, 10000) * 0.0 + 1.0
    plt.plot(Nd_all, omega_ratio_ref, '--k')
    plt.xlim(10**0, 2*10**5)
    plt.ylim(0.0, 1.02)
    plt.xscale("log")
    if save_fig:
        figname = "Freq_ratio_vs_decoys_Gamma_%1.1f_kappa_%1.1f_eps_%1.1f" % (Gamma, kappa, eps)
        figname = figname.replace(".", "_")
        figname = figname + ".pdf"
        save_dirname1 = "../Plots/"
        plt.savefig(save_dirname1 + figname)
        print("saved:", save_dirname1+figname)
    plt.show()

###############################################################################################
# Compute frequency ratio for  the repressilator
###############################################################################################
def freq_ratio(Np = 10.0, Nd = 100.0, beta = 10.0, eps = 10.0, kappa = 10.0):
    Np_r = Np*10**-2
    Nd_r = Nd*10**-2

    # kappa =  r_2 #alpha * sigma / (gamma_p * gamma_m)
    r_1 = 1.0
    r_2 = eps
    # solving equilibrium equation x + x^3 = kappa*N_p_r
    coeffs = [1.0, 0.0, 1.0, -kappa * Np_r]

    sols = np.roots(coeffs)
    bool_sols = np.isreal(sols)
    for i in range(3):
        if bool_sols[i]:
            des_sol = sols[i]
            break

    steady_x_real = des_sol.real
    p_x_steady = 1.0 + 4.0 * r_1 * steady_x_real + 4.0 * Np_r * steady_x_real / (
            1 + steady_x_real ** 2) ** 2 + 4.0 * Nd_r * r_2*steady_x_real / (
                         1 + r_2* steady_x_real ** 2) ** 2
    # solving equilibrium equation x + x^3 = kappa*N_p_r
    Np_r_ref = 1.0*10**-2
    coeffs0 = [1.0, 0.0, 1.0, -kappa * Np_r_ref]
    sols0 = np.roots(coeffs0)
    bool_sols0 = np.isreal(sols0)
    for i in range(3):
        if bool_sols0[i]:
            des_sol0 = sols0[i]
            break
    steady_x_real0 = des_sol0.real
    p_x_steady0 = 1.0 + 4.0 * r_1 * steady_x_real0 + 4.0 * Np_r * steady_x_real0 / (
            1 + steady_x_real0 ** 2) ** 2
    omega_r =  (beta+p_x_steady0)/(beta+p_x_steady)/Np_r*10**-2
    return omega_r

###############################################################################################
# Compute frequency of the repressilator
###############################################################################################
def repressilator_freq(Nd_r=0.0, beta = 10.0, eps = 1.0):

    r = 1.0
    Np_r = 1.0
    kappa = 500.0

    # solving equilibrium equation x + x^3 = kappa*N_p_r
    coeffs = [1.0, 0.0, 1.0, -kappa*Np_r]

    sols = np.roots(coeffs)
    bool_sols = np.isreal(sols)
    for i in range(3):
        if bool_sols[i]:
            des_sol = sols[i]
            break
    steady_x_real = des_sol.real
    p_x_steady = 1.0 + 4.0 * r * steady_x_real + 4.0 * Np_r * steady_x_real / (
                1 + steady_x_real ** 2) ** 2 + 4.0 * Nd_r*eps*steady_x_real / (
                             1 + eps* steady_x_real ** 2) ** 2

    omega = sqrt(3) * steady_x_real ** 3 / (kappa * Np_r) *beta / (beta + p_x_steady)
    return omega
###############################################################################################
#  Execute main()
###############################################################################################

if __name__ == "__main__":
    main()
