#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline

colors = mpl.cm.get_cmap('tab20').colors

temps = [91, 101, 111, 120, 130, 140, 150, 160, 170, 180, 190]

R = 8.3144598E-3
F0 = 0
Vi_op = 1.02E-25
Vi_cp = 4.736E-26

def get_nrtp(P, T, Nads):
    return Nads*((R*T)/P)

for n,temp in enumerate(temps):
    data = pd.read_pickle("./adsorption/summary_"+str(temp)+"K.pkl")


    nrtp_op = get_nrtp(data["pressure"], temp, data["amount_op"])
    nrtp_interp_op = InterpolatedUnivariateSpline(data["pressure"], nrtp_op, k=1)

    nrtp_cp = get_nrtp(data["pressure"], temp, data["amount_cp"])
    nrtp_interp_cp = InterpolatedUnivariateSpline(data["pressure"], nrtp_cp, k=1)

    pressure = np.geomspace(min(data["pressure"]), max(data["pressure"]), num=100000)
    
    os_op = np.array([-1*nrtp_interp_op.integral(0, p) for p in pressure])
    os_cp =  np.array([-1*nrtp_interp_cp.integral(0, p) for p in pressure])

    plt.semilogx(pressure,os_op-os_cp, color=colors[n], label=(str(temp)+"K"))
    # plt.plot(pressure,os_cp, '--', color=colors[n])
    


    # idx = np.argwhere(np.diff(np.sign(amount_cp_interp(pressure) - amount_op_interp(pressure)))).flatten()

    # print([temp, pressure[idx],amount_op_interp(pressure[idx])])
    # plt.semilogx(pressure, amount_op_interp(pressure), color=colors[n], alpha=0.3)
    # plt.plot(pressure, amount_cp_interp(pressure), '--', color=colors[n],  alpha=0.3)
    # plt.scatter(pressure[idx], amount_op_interp(pressure[idx]), color=colors[n], label=(str(temp)+"K"))


plt.ylim([-500,1200])
plt.ylabel("$\Delta f_i$ / kJ$\,$mol$^{-1}$")
plt.xlabel("pressure / kPa")
plt.legend(loc=3)
plt.savefig("./adsorption/osmotic_diff.png", bbox_inches="tight")
#plt.show()


