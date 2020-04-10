from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import socket
import time
import h5py
########################################################################################################################

# Parameters slightly different that those in MilJohWei2008 for toggle switch

# following parameters have the unit of min inverse
# alpha = 10.12  (assuming alpha_1 = alpha_2)
# beta = 0.22 (assuming beta_1 = beta_2)
# gamma_p = 0.05
# sigma = 0.5
# C = 10^-9 # mol inv
# gamma_m = 0.1

# following parameters are ratio (no units)
# cp = 10^7
# cd = 10^7

########################################################################################################################
# Main
#######################################################################################################################
def main():
    toggle_switch()

########################################################################################################################
# Basin of attraction
#######################################################################################################################

def toggle_switch():
    Np_all = np.array([1.0, 3.0])
    Nd_all = np.array([1000, 3000, 6000, 8000])
    eps = 1.0
    for Np in Np_all:
        for Nd in Nd_all:
            Nd1 = 100.0
            basin_attraction_identifier1(save_data=True, eps=eps, Np=Np, Nd1=Nd1, Nd2=Nd)
            Nd2 = 100.0
            basin_attraction_identifier1(save_data=True, eps=eps, Np=Np, Nd1=Nd, Nd2=Nd2)
########################################################################################################
# Identify the basin of attraction
########################################################################################################
def basin_attraction_identifier1(save_data=True, eps = 0.5, Np = 1.0, Nd1 = 100, Nd2 = 1000):
    if Np == 1:
        x1_all = np.linspace(0.01, 7.0, 500)
    else:
        x1_all = np.linspace(0.01, 18.0, 500)
    x1_manifold = []
    x2_manifold = []

    for x1_begin in x1_all:
        x1 = x1_begin
        x2_min_new = 0.01
        x2_max_new = x1_begin + 100.0
        for j in range(10):
            x2_min = x2_min_new
            x2_max = x2_max_new
            # from IPython import embed
            # embed()
            x2_max_new, x2_min_new, Error = find_search_interval(x1=x1,
                                                                 x2_min=x2_min, x2_max=x2_max ,eps = eps, Np = Np,
                                                                 Nd1 = Nd1, Nd2=Nd2)

            diff_x2 = abs(x2_min_new - x2_max_new)
            if Error==0:
                if diff_x2 < 10 **-7:
                    #print("found the point")
                    x2_manifold_point = x2_min_new
                    x1_manifold.append(x1)
                    x2_manifold.append(x2_manifold_point)
                    break
            elif Error == 1.0:
                x1_manifold.append(x1)
                x2_manifold.append(0.0)
                #print("appending 0.0 for x_2")
                break
        #print("--------")

    if save_data:
        dirname = "../Data/"
        fname = "Separatrices_Np_%d_Nd1_%d_Nd2_%d_eps%1.1f.h5" % (Np,Nd1,Nd2, eps)
        fname = fname.replace(".", "_")
        fname = fname + ".h5"
        hf = h5py.File(dirname + fname, 'w')
        hf.create_dataset('Gene1', data=x1_manifold)
        hf.create_dataset('Gene2', data=x2_manifold)
        hf.close()

#####################################################################################################
# Find search interval
#####################################################################################################
def find_search_interval(x1 = 1.0, x2_min = 1.0, x2_max = 2.0, Np=3.0, Nd1=100.0, Nd2 = 1000.0, eps=0.5):
    t = arange(0, 1000, 0.01)
    x2_all = np.linspace(x2_min,x2_max, 100)
    # from IPython import embed
    # embed()
    x2_min_for_return = x2_min
    x2_max_for_return = x2_max

    Error = 1.0
    check = 1.0
    for x2 in x2_all:
        vec = [x1, 0.0, x2, 0.0]
        #print(vec)
        sol = odeint(model_toggle_switch, vec, t, args=(Np, Nd1, Nd2, eps))
        diff_stead_val = sol[-1][0] - sol[-1][2]
        #print(x1,x2, diff_stead_val)
        if x2 > x2_min and check*diff_stead_val < 0:
            #print("change of sign detected")
            x2_max_for_return = x2
            Error = 0.0
            break
        else:
            check = diff_stead_val/abs(diff_stead_val)
            x2_min_for_return = x2
            # from IPython import embed
            # embed()

    if Error == 1.0:
        #print("point not found, returning zero...")
        x2_min_for_return = 0.0
        x2_max_for_return = 0.0
    return  x2_max_for_return, x2_min_for_return, Error


#####################################################################################################
# Model
#####################################################################################################
def model_toggle_switch(vec,t, Np, Nd1, Nd2, eps):

    x1 = vec[0]
    m1 = vec[1]
    x2 = vec[2]
    m2 = vec[3]

    # Parameters
    r1 = 1.0
    eps1 = eps
    eps2 = 1.0
    kappa1_1 = 6.325 * 100
    kappa1_2 = 0.1375 * 100
    kappa2_1 = 6.325 * 100
    kappa2_2 = 0.1375 * 100
    Gamma = 0.5
    Np_r = Np*10**-2
    Nd_r1 = Nd1*10**-2
    Nd_r2 = Nd2 * 10 ** -2
    p_x1 = 1.0 + 4.0 * r1 * x1 + 4.0 * Np_r * x1 / (1 + x1 ** 2) ** 2 + 4.0 * Nd_r1 * eps1 * x1 / (
                1 + eps1 * x1 ** 2) ** 2
    p_x2 = 1.0 + 4.0 * r1 * x2 + 4.0 * Np_r * x2 / (1 + x2 ** 2) ** 2 + 4.0 * Nd_r2 * eps2 * x2 / (
                1 + eps2 * x2 ** 2) ** 2

    dx1 = -Gamma * (x1 - m1)/p_x1
    dm1 = kappa1_1 * Np_r/(1 + x2 ** 2) + kappa1_2 * Np_r*x2**2/(1 + x2 ** 2) - m1

    dx2 = -Gamma * (x2 - m2)/p_x2
    dm2 = kappa2_1 * Np_r/(1 + x1 ** 2) + kappa2_2 * Np_r*x1**2/(1 + x1 ** 2) - m2

    return [dx1, dm1, dx2, dm2]


##############################################################################
# Main()
##############################################################################
if __name__ == "__main__":
    main()