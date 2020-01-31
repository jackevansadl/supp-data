#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from glob import glob
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.signal import argrelextrema

data_temp91 = pd.read_pickle("./stress/loadingstress_91K.pkl")
emptystress_91 = np.mean(data_temp91[data_temp91["nads"] == 0]["stress"])
data_temp91["relative_stress"] = (data_temp91["stress"]-emptystress_91)*0.101325
data_isotherm = pd.read_pickle("./adsorption/summary_91K.pkl")
amount_isotherm = data_isotherm["amount_op"]

f, axarr = plt.subplots(2, sharex=True)

average_stress = []
stddev_stress = []
for amount in amount_isotherm:
    df_closevalues = data_temp91.iloc[(data_temp91['nads']-amount).abs().argsort()[:5]]
    average_stress.append(np.mean(df_closevalues["relative_stress"]))
    stddev_stress.append(np.std(df_closevalues["relative_stress"]))
average_stress = np.array(average_stress)
stddev_stress = np.array(stddev_stress)

axarr[1].fill_between(data_isotherm["pressure"],average_stress-stddev_stress,average_stress+stddev_stress, color="C0", alpha=0.3)
axarr[1].semilogx(data_isotherm["pressure"],average_stress, 'o-', color="C0")

data_temp91 = pd.read_pickle("./stress/loadingstress_110K.pkl")
emptystress_91 = np.mean(data_temp91[data_temp91["nads"] == 0]["stress"])
data_temp91["relative_stress"] = (data_temp91["stress"]-emptystress_91)*0.101325
data_isotherm = pd.read_pickle("./adsorption/summary_111K.pkl")
amount_isotherm = data_isotherm["amount_op"]

average_stress = []
stddev_stress = []
for amount in amount_isotherm:
    df_closevalues = data_temp91.iloc[(data_temp91['nads']-amount).abs().argsort()[:5]]
    average_stress.append(np.mean(df_closevalues["relative_stress"]))
    stddev_stress.append(np.std(df_closevalues["relative_stress"]))
average_stress = np.array(average_stress)
stddev_stress = np.array(stddev_stress)

axarr[1].fill_between(data_isotherm["pressure"],average_stress-stddev_stress,average_stress+stddev_stress, color="C1", alpha=0.3)
axarr[1].semilogx(data_isotherm["pressure"],average_stress, 'o-', color="C1")

data_temp91 = pd.read_pickle("./stress/loadingstress_130K.pkl")
emptystress_91 = np.mean(data_temp91[data_temp91["nads"] == 0]["stress"])
data_temp91["relative_stress"] = (data_temp91["stress"]-emptystress_91)*0.101325
data_isotherm = pd.read_pickle("./adsorption/summary_130K.pkl")
amount_isotherm = data_isotherm["amount_op"]

average_stress = []
stddev_stress = []
for amount in amount_isotherm:
    df_closevalues = data_temp91.iloc[(data_temp91['nads']-amount).abs().argsort()[:5]]
    average_stress.append(np.mean(df_closevalues["relative_stress"]))
    stddev_stress.append(np.std(df_closevalues["relative_stress"]))
average_stress = np.array(average_stress)
stddev_stress = np.array(stddev_stress)

axarr[1].fill_between(data_isotherm["pressure"],average_stress-stddev_stress,average_stress+stddev_stress, color="C2", alpha=0.3)
axarr[1].semilogx(data_isotherm["pressure"],average_stress, 'o-', color="C2")



axarr[1].set_ylim([-40,40])

temps = [91, 111, 130]

for n,temp in enumerate(temps):
    data = pd.read_pickle("./adsorption/summary_"+str(temp)+"K.pkl")

    amount_op_interp = interp1d(data["pressure"], data["amount_op"], fill_value="extrapolate")
    amount_cp_interp = interp1d(data["pressure"], data["amount_cp"], fill_value="extrapolate")

    pressure = np.geomspace(min(data["pressure"]), max(data["pressure"]), num=100000)

    idx = np.argwhere(np.diff(np.sign(amount_cp_interp(pressure) - amount_op_interp(pressure)))).flatten()

    print([temp, pressure[idx],amount_op_interp(pressure[idx])])
    axarr[0].plot(pressure, amount_op_interp(pressure), color=("C"+str(n)), alpha=0.3)
    axarr[0].plot(pressure, amount_cp_interp(pressure), '--', color=("C"+str(n)),  alpha=0.3)
    axarr[0].scatter(pressure[idx], amount_op_interp(pressure[idx]), color=("C"+str(n)), label=(str(temp)+"K"))
    axarr[1].axvline(x=pressure[idx], color=("C"+str(n)))
    print(amount_op_interp(pressure[idx]))






axarr[1].set_ylabel("stress / MPa")
# axarr[0].set_ylabel("amount adsorbed / molec. per unit cell")
axarr[1].set_xlabel("pressure / kPa")
plt.savefig("./stress/stress_isotherm.png", bbox_inches="tight")


#%%
