from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import arange
# import CPG
import socket
import time

###############################################################################################
# Main
###############################################################################################
def main():
    """
    Interesting parameters:


    :return:
    """
    # Regime C # check the value of Np
    beta = 0.8
    kappa = 600
    eps = 1.0
    repressilator(beta=beta,kappa=kappa, eps=eps, Np = 3.0, save_fig=True)

    # Regime A
    beta =  12.0
    kappa = 145.0
    eps = 1
    repressilator(beta=beta,kappa=kappa, eps=eps, Np = 3.0, save_fig=True)

    # Increased amplitude # put Np = 1
    beta = 40
    kappa = 800
    eps = 1
    repressilator(beta=beta, kappa=kappa, eps=eps, Np = 1.0, save_fig=True)


###############################################################################################
# Time evolution
###############################################################################################
def repressilator(beta = 10.0, kappa = 5.0, eps = 0.1, Np = 1, save_fig=False):
    t = arange(0, 5000, 0.01)
    vec = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    if beta == 40:
        Nd_all = np.array([3000.0])
    else:
        Nd_all = np.array([1000.0])

    for Nd in Nd_all:
        Nd1 = 0.0
        sol1 = odeint(model_repressilator, vec, t, args=(Np,Nd1, beta, kappa, eps))
        sol2 = odeint(model_repressilator, vec, t, args=(Np,Nd, beta, kappa, eps))
        plt.figure()
        plt.plot(t, sol2[:, 0], 'r')
        plt.plot(t, sol1[:, 0])
        plt.ylabel("Monomer concentration")
        plt.xlabel("rescaled time")
        if beta == 0.8:
            figname = r"Monomer_conc_repressilator_Gamma_%1.1f_kappa_%1.1f_Np_%d_Nd_%d" % (beta, kappa, Np, Nd)
            plt.xlim([0, 3000])
        elif beta == 12:
            figname = r"Monomer_conc_repressilator_Gamma_%1.1f_kappa_%1.1f_Np_%d_Nd_%d" % (beta, kappa, Np, Nd)
            plt.xlim([0, 400])
        elif beta==40:
            figname = r"Monomer_conc_repressilator_Gamma_%1.1f_kappa_%1.1f_Np_%d_Nd_%d" % (beta, kappa, Np, Nd)
            plt.xlim([0, 100])
        if save_fig:
            figname = figname.replace(".", "_")
            figname = figname + ".pdf"
            save_dirname1 = "../Plots/Large_decoy_amp/"
            plt.savefig(save_dirname1 + figname)
            plt.show()
            print("saved:", figname)
        plt.show()



###############################################################################################
# Model_repressilator
###############################################################################################

def model_repressilator(vec,t, Np, Nd, beta, kappa, eps):

    x1 = vec[0]
    m1 = vec[1]
    x2 = vec[2]
    m2 = vec[3]
    x3 = vec[4]
    m3 = vec[5]

    # Parameters
    r1 = 1.0
    r2 = eps

    Np_r = Np * 10**-2
    Nd_r = Nd * 10**-2

    p_x1 = 1.0 + 4.0*r1*x1 + 4.0*Np_r*x1/(1+x1**2)**2 + 4.0*Nd_r*r2*x1/(1+r2*x1**2)**2
    p_x2 = 1.0 + 4.0*r1*x2 + 4.0*Np_r*x2/(1+x2**2)**2 + 4.0*Nd_r*r2*x2/(1+r2*x2**2)**2
    p_x3 = 1.0 + 4.0*r1*x3 + 4.0*Np_r*x3/(1+x3**2)**2 + 4.0*Nd_r*r2*x3/(1+r2*x3**2)**2


    # from IPython import embed
    # embed()

    # p_x = 1 + 4.0*cp + 4.0*Np*(cp*cd1*x)/(1+cp*cd1*x**2)**2 + 4.0*Nd*(cp*x)/(cp*x**2)**2 # k- = 0 (perfect protection)

    dx1 = -beta * (x1 - m1)/p_x1
    dm1 = kappa * Np_r/(1 + x3 ** 2) - m1

    dx2 = -beta * (x2 - m2)/p_x2
    dm2 = kappa * Np_r/(1 + x1 ** 2) - m2

    dx3 = -beta * (x3 - m3)/p_x3
    dm3 = kappa * Np_r/(1 + x2 ** 2) - m3

    return [dx1, dm1, dx2, dm2, dx3, dm3]

###############################################################################################
# Send Email
###############################################################################################
def send_email(script_name):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    from email.mime.base import MIMEBase
    # creates SMTP session
    s = smtplib.SMTP('smtp.gmail.com', 587)

    # start TLS for security
    s.starttls()

    # Authentication
    s.login("sdas.script.info@gmail.com", "pythonscript2019")

    # message to be sent
    msg = MIMEMultipart()
    # msg['From'] = 'sdas.script.info@gmail.com'
    # msg['To'] = 'swetdas@gmail.com'
    msg['Subject'] = 'Amplitude computations finished'

    body = script_name
    msg.attach(MIMEText(body, 'plain'))

    message = msg.as_string()


    # sending the mail
    s.sendmail("sdas.script.info@gmail.com", "swetdas@gmail.com", message)

    # terminating the session
    s.quit()


###############################################################################################
###############################################################################################

if __name__ == "__main__":
    main()