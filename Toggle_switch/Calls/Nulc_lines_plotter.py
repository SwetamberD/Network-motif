from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
# import CPG

########################################################################################################################
def main():
########################################################################################################################
    Np_all = np.array([1.0, 3.0])
    # Plot nulc lines

    for Np in Np_all:
        sol1_x, sol1_y, sol2_x, sol2_y = nulc_lines(Np = Np)
        plt.plot(sol1_x, sol1_y, 'grey')
        plt.plot(sol2_x, sol2_y, 'grey')

        # Plot separatix in the absence of decoy sites
        if Np == 1.0:
            x_sep = np.linspace(0, 7, 10)
            y_sep = x_sep
            plt.plot(x_sep, y_sep, '--k')
            plt.xlim([0.0, 7.0])
            plt.ylim([0.0, 7.0])
        else:
            x_sep = np.linspace(0, 17, 10)
            y_sep = x_sep
            plt.plot(x_sep, y_sep, '--k')
            plt.xlim([0.0, 16.50])
            plt.ylim([0.0, 16.50])
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

        figname = "Nulc_lines_for_Np_%d.pdf" % Np
        plt.savefig("../Plots/"+figname)
        print("saved:",figname)
        plt.show()

########################################################################################################################
def nulc_lines(Np=1.0):
########################################################################################################################
    k1 = 6.325 * 100
    k2 = 0.1375 * 100
    Np = Np* 10**-2 # due to scaling
    x2_all = np.linspace(0.1, 20, 1000)
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
