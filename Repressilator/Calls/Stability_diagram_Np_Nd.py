from __future__ import division, print_function

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import multiprocessing
import time
#import CPG
from sympy import *
from sympy import Matrix

##############################################################################################
# main
###############################################################################################

def main():
    comparison_plot_examples(save_fig=True)

##############################################################################################
# Compare stability diagram
###############################################################################################

def comparison_plot(Np = 1.0,Nd1 = 0.0, Nd2 = 10.0, save_fig=False):
    eps = 1.0
    Np0, beta1_0, beta2_0 = stability_diagram(Np_val = Np, Nd_val=Nd1, eps= eps)
    Np1, beta1_1, beta2_1 = stability_diagram(Np_val = Np, Nd_val=Nd2, eps=eps)
    plt.figure()
    plt.plot(Np0, beta1_0, 'b')
    plt.plot(Np0, beta2_0, 'b')
    plt.plot(Np1, beta1_1, 'r')
    plt.plot(Np1, beta2_1, 'r')
    plt.fill_between(Np0, beta2_0, beta1_1, where=beta2_0>beta1_1, color="grey", alpha=0.5)  # hatch = '|' )
    plt.fill_between(Np1, beta2_0, beta1_1, where=beta2_0<beta1_1, color='white', alpha=0.5)
    plt.xlabel(r"$\kappa$")
    plt.ylabel(r"$\Gamma$")
    plt.rcParams.update({'font.size': 10})
    plt.xlim([6*10, 10**3])
    #plt.ylim([1,100])
    plt.yscale("log")
    plt.xscale("log")
    title = r"$N_p$ = %d, $N_d$ = %d" %(Np, Nd2)
    plt.title(title,fontdict = {'fontsize' : 8})
    if save_fig:
        figname = "Stability_diagram_comparison_Np_1_and_3_Nd_%d" % Nd
        figname = figname.replace(".","_")
        figname = figname + ".pdf"
        save_dirname1 = "../Plots/"
        plt.savefig(save_dirname1 + figname)
    plt.show()


##############################################################################################
# Compare stability diagram (three examples)
###############################################################################################

def comparison_plot_examples(save_fig=False):
    eps = 0.1
    Np_val1 = 1.0
    Np_val2 = 3.0
    Nd1 = 0.0
    Nd2 = 1000
    kappa0, beta1_0, beta2_0 = stability_diagram(Np_val = Np_val1, Nd_val=Nd1, eps= eps)
    kappa1, beta1_1, beta2_1 = stability_diagram(Np_val = Np_val1, Nd_val=Nd2, eps=eps)
    kappa2, beta1_2, beta2_2 = stability_diagram(Np_val = Np_val2, Nd_val=Nd1, eps=eps)
    kappa3, beta1_3, beta2_3 = stability_diagram(Np_val = Np_val2, Nd_val=Nd2, eps=eps)
    plt.figure()
    plt.plot(kappa0, beta1_0, 'g')
    plt.plot(kappa0, beta2_0, 'g')
    plt.plot(kappa1, beta1_1, 'r')
    plt.plot(kappa1, beta2_1, 'r')
    plt.plot(kappa2, beta1_2, 'g')
    plt.plot(kappa2, beta2_2, 'g')
    plt.plot(kappa2, beta1_3, 'r')
    plt.plot(kappa3, beta2_3, 'r')
    # plt.fill_between(Np0, beta2_0, beta1_1, where=beta2_0>beta1_1, color="grey", alpha=0.5)  # hatch = '|' )
    # plt.fill_between(Np1, beta2_0, beta1_1, where=beta2_0<beta1_1, color='white', alpha=0.5)
    plt.xlabel(r"$\kappa$")
    plt.ylabel(r"$\Gamma$")
    plt.rcParams.update({'font.size': 10})
    #plt.text(10, 300,r"$A$")
    # plt.text(8, 100,r"$A$")
    # plt.text(10, 2, r"$C$")
    # plt.text(100, 20, r"$B$")
    #plt.xlim([6*10, 10**3])
    #plt.ylim([1,100])
    plt.yscale("log")
    plt.xscale("log")
    # title = r"$N_p$ = %d, $N_d$ = %d" %(Np, Nd2)
    # plt.title(title,fontdict = {'fontsize' : 8})
    if save_fig:
        figname = "Stability_diagram_comparison_Np_Nd_1_1000_3_1000"
        figname = figname.replace(".","_")
        figname = figname + ".pdf"
        save_dirname1 = "../Plots/"
        #save_dirname2 = "../../../Fig_files/Fig_Repressilator/"
        plt.savefig(save_dirname1 + figname)
        #plt.savefig(save_dirname2 + figname)
        print("saved:", figname)
    else:
        plt.show()


##############################################################################################
# Stability diagram
###############################################################################################
def stability_diagram(Np_val = 1.0, Nd_val = 0.0, eps = 10.1):
    beta1_all = []
    beta2_all = []
    kappa_sol_all = []
    Np = Np_val * 10**-2

    kappa_all = np.linspace(1, 1000, 100000)
    Nd_r = Nd_val * 10**-2

    for kappa in kappa_all:

        coeffs = [1.0, 0.0, 1.0, -kappa * Np]

        sols = np.roots(coeffs)
        bool_sols = np.isreal(sols)
        for i in range(3):
            if bool_sols[i]:
                des_sol = sols[i]
                break
        steady_x_real = des_sol.real
        A = - 2.0*steady_x_real**3/(kappa*Np)
        #print(A)
        beta1 = (3.0*A**2 - 4.0*A - 8.0)/(4.0*A + 8.0) + A*np.sqrt(9.0*A**2 - 24.0*A - 48.0)/(4.0*A + 8.0)
        beta2 = (3.0*A**2 - 4.0*A - 8.0)/(4.0*A + 8.0) - A*np.sqrt(9.0*A**2 - 24.0*A - 48.0)/(4.0*A + 8.0)

        r_1 = 1.0
        r_2 = eps
        Np_r = Np

        p_x_steady = 1.0 + 4.0 * r_1 * steady_x_real + 4.0 * Np_r * steady_x_real / (
                1 + steady_x_real ** 2) ** 2 + 4.0 * Nd_r * r_2*steady_x_real / (
                             1 + r_2 * steady_x_real ** 2) ** 2

        #print(beta1)

        beta1_p = beta1 *p_x_steady
        beta2_p = beta2 *p_x_steady
        if Np == 100:
            print(1000.0/p_x_steady)

        if beta1.imag == 0:
            beta1_all.append(beta1_p)
            beta2_all.append(beta2_p)
            kappa_sol_all.append(kappa)
    return np.array(kappa_sol_all), np.array(beta1_all), np.array(beta2_all)


##############################################################################################
###############################################################################################

if __name__ == "__main__":
    main()
