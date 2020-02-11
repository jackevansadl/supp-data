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
data_temp110 = pd.read_pickle("./stress/loadingstress_110K.pkl")
data_temp130 = pd.read_pickle("./stress/loadingstress_130K.pkl")

emptystress_91 = np.mean(data_temp91[data_temp91["nads"] == 0]["stress"])
emptystress_110 = np.mean(data_temp110[data_temp110["nads"] == 0]["stress"])
emptystress_130 = np.mean(data_temp130[data_temp130["nads"] == 0]["stress"])

data_temp91["nads_nga"] = data_temp91["nads"]-484
data_temp91["relative_stress"] = (data_temp91["stress"]-emptystress_91)*0.101325

data_temp110["nads_nga"] = data_temp110["nads"]-461
data_temp110["relative_stress"] = (data_temp110["stress"]-emptystress_110)*0.101325

data_temp130["nads_nga"] = data_temp130["nads"]-430
data_temp130["relative_stress"] = (data_temp130["stress"]-emptystress_130)*0.101325

data_temp91 = data_temp91[data_temp91['nads_nga'].between(-100, 300, inclusive=True)]
data_temp110 = data_temp110[data_temp110['nads_nga'].between(-100, 300, inclusive=True)]
data_temp130 = data_temp130[data_temp130['nads_nga'].between(-100, 300, inclusive=True)]

data_temp91 = data_temp91.sort_values(by=['nads_nga'])
data_temp110 = data_temp110.sort_values(by=['nads_nga'])
data_temp130 = data_temp130.sort_values(by=['nads_nga'])


# plt.plot(data_temp91["nads"]-484,(data_temp91["stress"]-emptystress_91)*0.101325, 'o', label="91K")
# plt.plot(data_temp110["nads"]-461,(data_temp110["stress"]-emptystress_110)*0.101325, 'o', label="110K")
# plt.plot(data_temp130["nads"]-430,(data_temp130["stress"]-emptystress_130)*0.101325, 'o', label="130K")


plt.plot(data_temp91["nads_nga"], data_temp91["relative_stress"], 'o', label="91K")
plt.plot(data_temp110["nads_nga"], data_temp110["relative_stress"], 'o', label="110K")
plt.plot(data_temp130["nads_nga"], data_temp130["relative_stress"], 'o', label="130K")

xp = np.linspace(-400, 1000, 100)
p2 = np.poly1d(np.polyfit(data_temp91["nads_nga"], data_temp91["relative_stress"], 2))
plt.plot(xp, p2(xp), color='C0')
p2 = np.poly1d(np.polyfit(data_temp110["nads_nga"], data_temp110["relative_stress"], 2))
plt.plot(xp, p2(xp), color='C1')
p2 = np.poly1d(np.polyfit(data_temp130["nads_nga"], data_temp130["relative_stress"], 2))
plt.plot(xp, p2(xp), color='C2')



plt.axhline(y=-30, color="K")
# plt.axvline(x=461, color="C1")
# plt.axvline(x=430, color="C2")
plt.ylim([-60, 0])
plt.xlim([-50, 300])
plt.legend()
plt.xlabel("N$_{\mathrm{NGA}}$ / molec. per unit cell")
plt.ylabel("stress / MPa")
plt.savefig("./stress/stress_loading_NGA.png", bbox_inches="tight")

#%%

