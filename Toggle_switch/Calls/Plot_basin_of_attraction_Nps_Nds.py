from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
import h5py

########################################################################################################################
# main
########################################################################################################################
def main():
    save_fig = True
    open_dir = "../Data/"
    Np_all = np.array([1.0, 3.0])
    Nd_all = np.array([3000, 8000])
    inds = np.array([1.0,2.0])
    eps = 1.0
    for Np in Np_all:
        for Nd in Nd_all:
            for ind in inds:
                if ind == 1:
                    Nd1 = Nd
                    Nd2 = 100
                    fname = "Separatrices_Np_%d_Nd1_%d_Nd2_%d_eps%1.1f.h5" % (Np,Nd1,Nd2, eps)
                    alpha = 1.0
                else:
                    Nd1 = Nd
                    Nd2 = 100
                    fname = "Separatrices_Np_%d_Nd1_%d_Nd2_%d_eps%1.1f.h5" % (Np, Nd1, Nd2, eps)
                    alpha = 0.3
                fname = fname.replace(".", "_")
                fname = fname + ".h5"

                hf = h5py.File(open_dir + fname, 'r')
                Gene1= hf.get('Gene1')
                Gene1 = np.array(Gene1)
                Gene2 = hf.get('Gene2')
                Gene2 = np.array(Gene2)
                Half_times = hf.get('Half_times')
                Half_times = np.array(Half_times)
                hf.close()
                label = r"$N_p$ = %d, $N_{d1}$ = %d, $N_{d2}$ = %d" % (Np, Nd1, Nd2)
                if Nd == 3000:
                    plt.plot(Gene1, Gene2, 'm', label=label, alpha = alpha)
                else:
                    plt.plot(Gene1, Gene2, 'r', label=label, alpha = alpha)


        # Plot nulc lines
        sol1_x, sol1_y, sol2_x, sol2_y = nulc_lines(Np=Np)
        plt.plot(sol1_x, sol1_y, 'grey')
        plt.plot(sol2_x, sol2_y, 'grey')

        # Plot separatix in the absence of decoy sites
        x_sep = np.linspace(0, 20, 10)
        y_sep = x_sep
        plt.plot(x_sep,y_sep,'--k')
        plt.xlim([0.0, 7.0])
        plt.ylim([0.0, 7.0])
        plt.xlabel("Gene 1")

        plt.ylabel("Gene 2")

        if Np == 1:
            plt.xlim([0.0, 7.0])
            plt.ylim([0.0, 7.0])
            text_1 = r"(Low, High)"
            text_2 = r"(High, Low)"
            plt.text(5.6, 0.6, text_1)
            plt.text(0.5, 5.6, text_2)

        else:
            plt.xlim([0.0, 18.0])
            plt.ylim([0.0, 18.0])
            text_1 = r"(Low, High)"
            text_2 = r"(High, Low)"
            plt.text(15.0, 2.0, text_1)
            plt.text(2.0, 15.0, text_2)

        if save_fig == True:
            figname = "Basin_attraction_for_different_Nds_for_Np_%d.pdf" % Np
            save_dir = "../Plots/"
            plt.savefig(save_dir + figname)
            print("saved:", figname)
        plt.show()
########################################################################################################################
# Nulc lines
########################################################################################################################
def nulc_lines(Np=1):
    k1 = 6.325 * 100
    k2 = 0.1375 * 100


    Np = Np* 10**-2
    x2_all = np.linspace(0.1, 25, 1000)

    sol1_x =[]
    sol1_y =[]
    sol2_x =[]
    sol2_y =[]

    for x2 in x2_all:
        f1 = k1*Np/(1+x2**2) + k2*Np*x2**2/(1+x2**2)
        f2 = np.sqrt((k1*Np- x2)/(x2 - k2*Np))

        if f2.imag == 0:
            sol1_y.append(f1)
            sol1_x.append(x2)
            sol2_y.append(f2)
            sol2_x.append(x2)

    for i in range(len(sol1_x)):
        error = abs(sol1_y[i] - sol2_y[i])
        if error < 10**-3:
            print(sol1_x[i], sol1_y[i], sol1_y[i])

    return sol1_x, sol1_y, sol2_x, sol2_y

########################################################################################################################

if __name__ == "__main__":
    main()
