from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import h5py
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
    half_time_computation(save_data = True)

########################################################################################################################
def half_time_computation(save_data=True):
    t_max = 8000
    t_step = 0.01
    n_times = int(t_max / t_step)
    t = arange(0, t_max, t_step)
    Np_all = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    Nd_all = np.array([0.0])
    Half_life_for_Np = []
    for Np in Np_all:
        # find the stable fixed point
        fixed_point = stable_fixed_point(Np=Np)
        # generate 100 random initial conditions around the fixed point
        low_bound = fixed_point - 10**-2 # neighborhood of the fixed point
        up_bound = fixed_point + 10**-2 # neighborhood of the fixed point
        in_cods = np.random.uniform(low_bound, up_bound,100)
        avg_half_times = []
        for Nd in Nd_all:
            half_times = []
            for ic in in_cods:
                vec = [ic, ic]
                sol= odeint(model_positive_feedback, vec, t, args=(Np, Nd,))
                x_vals = sol[:, 0]
                steady_val = x_vals[-1]
                steady_state_check = abs(fixed_point - steady_val)
                if steady_state_check < 10 ** -3:
                    # determine the mid value
                    if x_vals[0] > fixed_point:
                        mid_stady_state = fixed_point + abs(x_vals[-1] - x_vals[0])/2.0
                    else:
                        mid_stady_state = fixed_point - abs(x_vals[-1] - x_vals[0])/2.0
                    diff_x_vals = x_vals - mid_stady_state
                    for i in range(n_times-1):
                        if diff_x_vals[i] *diff_x_vals[i+1] < 0:
                            # from IPython import embed
                            # embed()
                            half_times.append(i)
                            #print(i*0.01)
                            break
            half_time = np.mean(np.array(half_times) * 0.01)
            avg_half_times.append(half_time)
            Half_life_for_Np.append(half_time)
    Half_life_for_Np = np.array(Half_life_for_Np)

    if save_data:
        dirname = "../Data/"
        fname = "Half_times_for_different_Nps_Nd_0.h5"
        hf = h5py.File(dirname+fname, 'w')
        hf.create_dataset('Np', data=Np_all)
        hf.create_dataset('Half_times', data=Half_life_for_Np)
        hf.close()
        print("saved:", fname)

########################################################################################################################
# Steady state finder
########################################################################################################################
def stable_fixed_point(Np = 1.0):
    # Parameters
    kappa1 = 0.25*100
    kappa2 = 1.3*100
    Np_r =  Np* 10**-2 # due to rescaling
    coeffs = [1.0, -kappa2*Np_r, 1.0, -kappa1*Np_r]
    sols = np.roots(coeffs)
    bool_sols = np.isreal(sols)
    # print(sols)
    for i in range(3):
        if bool_sols[i]:
            des_sol = sols[i]
            break
    steady_x_real = des_sol.real
    return steady_x_real
########################################################################################################################
# Model
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

    Np_r = Np*0.01 # due to rescaling
    Nd_r = Nd*0.01 # due to rescaling


    p_x  = 1.0 + 4.0*r1*x + 4.0*Np_r*x/(1+x**2)**2 + 4.0*Nd_r*r2*x/(1+r2*x**2)**2

    dx = -beta * (x - m)/p_x
    dm = kappa1 * Np_r/(1 + x** 2) + kappa2 * Np_r*x** 2/(1 + x** 2)  - m

    return [dx, dm]

########################################################################################################################
if __name__ == "__main__":
    main()