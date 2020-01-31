#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from scipy.interpolate import interp1d

colors = mpl.cm.get_cmap('tab20').colors

temps = [91, 101, 111, 120, 130, 140, 150, 160, 170, 180, 190]

for n,temp in enumerate(temps):
    data = pd.read_pickle("./adsorption/summary_"+str(temp)+"K.pkl")

    amount_op_interp = interp1d(data["pressure"], data["amount_op"], fill_value="extrapolate")
    amount_cp_interp = interp1d(data["pressure"], data["amount_cp"], fill_value="extrapolate")

    pressure = np.geomspace(min(data["pressure"]), max(data["pressure"]), num=100000)

    idx = np.argwhere(np.diff(np.sign(amount_cp_interp(pressure) - amount_op_interp(pressure)))).flatten()

    print([temp, pressure[idx],amount_op_interp(pressure[idx])])
    plt.semilogx(pressure, amount_op_interp(pressure), color=colors[n], alpha=0.3)
    plt.plot(pressure, amount_cp_interp(pressure), '--', color=colors[n],  alpha=0.3)
    plt.scatter(pressure[idx], amount_op_interp(pressure[idx]), color=colors[n], label=(str(temp)+"K"))

plt.ylabel("amount adsorbed / molec. per unit cell")
plt.xlabel("pressure / kPa")
plt.legend()
#plt.savefig("./adsorption/crossingpoint.png", bbox_inches="tight")
plt.show()


