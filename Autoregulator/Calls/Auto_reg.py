from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
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
# Main
########################################################################################################################
def main():
    pos_feedback_change_Nd(save_fig=True)


########################################################################################################################
def pos_feedback_change_Nd(save_fig=False):
    t_max = 7000
    t_step = 0.01
    t = arange(0, t_max, t_step)
    vec = [0.1, 0.0]
    Nd_all = np.array([0.0, 500.0, 2000.0])
    Np = 1.0
    fixed_point = stable_fixed_point(Np=Np)
    print(r"steady state at:", fixed_point)
    for Nd in Nd_all:
        sol = odeint(model_positive_feedback, vec, t, args=(Np, Nd,))
        x_vals = sol[:, 0]
        plt.plot(t, x_vals, label= r"$N_d =$ %d"% Nd)
    plt.legend(loc='lower right')
    plt.xlabel('Rescaled time')
    plt.xlim([0, 700])
    plt.ylabel('Rescaled monomer concentration')
    if save_fig:
        figname = "monomer_conc_with_decoys_Np_%d.pdf" % Np
        dirname = "../Plots/"
        plt.savefig(dirname + figname)
        print("saved:",figname)
    plt.show()

########################################################################################################################
# Steady state finder
########################################################################################################################
def stable_fixed_point(Np=1.0):
    # Parameters
    kappa1 = 0.25*100
    kappa2 = 1.3*100
    Np_r =  Np* 10**-2 # due to scaling
    coeffs = [1.0, -kappa2*Np_r, 1.0, -kappa1*Np_r]
    sols = np.roots(coeffs)
    bool_sols = np.isreal(sols)
    for i in range(3):
        if bool_sols[i]:
            des_sol = sols[i]
            break
    steady_x_real = des_sol.real
    return steady_x_real

########################################################################################################################
def model_positive_feedback(vec, t, Np, Nd):

    x = vec[0]
    m = vec[1]

    # Parameters
    kappa1 = 0.25*100
    kappa2 = 1.3*100

    # # Parameters
    # kappa1 = 20.25 * 100
    # kappa2 = 10.7 * 100

    beta = 0.5

    r1 = 1.0
    r2 = 1.0 #eps

    Np_r = Np*.01
    Nd_r = Nd*.01


    p_x  = 1.0 + 4.0*r1*x + 4.0*Np_r*x/(1+x**2)**2 + 4.0*Nd_r*r2*x/(1+r2*x**2)**2

    dx = -beta * (x - m)/p_x
    dm = kappa1 * Np_r/(1 + x** 2) + kappa2 * Np_r*x** 2/(1 + x** 2)  - m

    return [dx, dm]

########################################################################################################################

if __name__ == "__main__":
    main()