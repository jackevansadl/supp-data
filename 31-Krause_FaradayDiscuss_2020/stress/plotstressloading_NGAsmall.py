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

mpl.rcParams.update({'font.size':10})
mpl.rcParams['figure.figsize'] = 4,3

data_temp91_big = pd.read_pickle("./stress/loadingstress_91K.pkl")
data_temp110_big = pd.read_pickle("./stress/loadingstress_110K.pkl")
data_temp130_big = pd.read_pickle("./stress/loadingstress_130K.pkl")

data_temp91 = pd.read_pickle("./stress/loadingstress_small_91K.pkl")
data_temp110 = pd.read_pickle("./stress/loadingstress_small_111K.pkl")
data_temp130 = pd.read_pickle("./stress/loadingstress_small_130K.pkl")

emptystress_91 = np.mean(data_temp91_big[data_temp91_big["nads"] == 0]["stress"])
emptystress_110 = np.mean(data_temp110_big[data_temp110_big["nads"] == 0]["stress"])
emptystress_130 = np.mean(data_temp130_big[data_temp130_big["nads"] == 0]["stress"])



data_temp91["nads_nga"] = data_temp91["nads"]-484
data_temp91["relative_stress"] = (data_temp91["stress"]-emptystress_91)*0.101325

data_temp110["nads_nga"] = data_temp110["nads"]-461
data_temp110["relative_stress"] = (data_temp110["stress"]-emptystress_110)*0.101325

data_temp130["nads_nga"] = data_temp130["nads"]-430
data_temp130["relative_stress"] = (data_temp130["stress"]-emptystress_130)*0.101325

# data_temp91 = data_temp91[data_temp91['nads_nga'].between(-100, 300, inclusive=True)]
# data_temp110 = data_temp110[data_temp110['nads_nga'].between(-100, 300, inclusive=True)]
# data_temp130 = data_temp130[data_temp130['nads_nga'].between(-100, 300, inclusive=True)]

data_temp91 = data_temp91.sort_values(by=['nads_nga'])
data_temp110 = data_temp110.sort_values(by=['nads_nga'])
data_temp130 = data_temp130.sort_values(by=['nads_nga'])

data_temp91.to_excel("./stress/loadingstress_small_91K.xlsx")
data_temp110.to_excel("./stress/loadingstress_small_110K.xlsx")
data_temp130.to_excel("./stress/loadingstress_small_130K.xlsx")

# plt.plot(data_temp91_big["nads"],(data_temp91_big["stress"]-emptystress_91)*0.101325, 'o', label="91K")
# plt.plot(data_temp110_big["nads"],(data_temp110_big["stress"]-emptystress_110)*0.101325, 'o', label="110K")
# plt.plot(data_temp130_big["nads"],(data_temp130_big["stress"]-emptystress_130)*0.101325, 'o', label="130K")

# plt.plot(data_temp91["nads_nga"], data_temp91["stress"], 'o', label="91K", color="C0", alpha=0.5)
# plt.plot(data_temp110["nads_nga"], data_temp110["stress"], 'o', label="110K", color="C1", alpha=0.5)
# plt.plot(data_temp130["nads_nga"], data_temp130["stress"], 'o', label="130K", color="C2", alpha=0.5)

uniquenads_91 = []
meanstress_91 = []
stdstress_91 = []
for nads in data_temp91["nads_nga"].unique():
    uniquenads_91.append(nads)
    data_unique = data_temp91[data_temp91["nads_nga"] == nads]
    meanstress_91.append(data_unique["relative_stress"].mean())
    stdstress_91.append(data_unique["relative_stress"].std())


uniquenads_111 = []
meanstress_111 = []
stdstress_111 = []
for nads in data_temp110["nads_nga"].unique():
    uniquenads_111.append(nads)
    data_unique = data_temp110[data_temp110["nads_nga"] == nads]
    meanstress_111.append(data_unique["relative_stress"].mean())
    stdstress_111.append(data_unique["relative_stress"].std())

uniquenads_130 = []
meanstress_130 = []
stdstress_130 = []
for nads in data_temp130["nads_nga"].unique():
    uniquenads_130.append(nads)
    data_unique = data_temp130[data_temp130["nads_nga"] == nads]
    meanstress_130.append(data_unique["relative_stress"].mean())
    stdstress_130.append(data_unique["relative_stress"].std())

plt.fill_between(uniquenads_91,meanstress_91-np.array(stdstress_91),meanstress_91+np.array(stdstress_91),color="C0",alpha=0.5)
plt.fill_between(uniquenads_111,meanstress_111-np.array(stdstress_111),meanstress_111+np.array(stdstress_111),color="C1",alpha=0.5)
plt.fill_between(uniquenads_130,meanstress_130-np.array(stdstress_130),meanstress_130+np.array(stdstress_130),color="C2",alpha=0.5)
plt.plot(uniquenads_91,meanstress_91,'-',color="C0", label="91 K")
plt.plot(uniquenads_111,meanstress_111,'-', color="C1", label="111 K")
plt.plot(uniquenads_130,meanstress_130,'-', color="C2", label="130 K")






# xp = np.linspace(-400, 1000, 100)
# p2 = np.poly1d(np.polyfit(data_temp91["nads_nga"], data_temp91["relative_stress"], 2))
# plt.plot(xp, p2(xp), color='C0')
# p2 = np.poly1d(np.polyfit(data_temp110["nads_nga"], data_temp110["relative_stress"], 2))
# plt.plot(xp, p2(xp), color='C1')
# p2 = np.poly1d(np.polyfit(data_temp130["nads_nga"], data_temp130["relative_stress"], 2))
# plt.plot(xp, p2(xp), color='C2')


plt.axhline(y=-20.5, color="K")
plt.axhspan(-20.5, -50, color="K", alpha=0.1)
# plt.axvline(x=461, color="C1")
# plt.axvline(x=430, color="C2")
# plt.xlim([0, 1200])
plt.ylim([-33, -10])
plt.legend()
plt.xlabel("N$_{\mathrm{NGA}}$ / molec. per unit cell")
plt.ylabel("stress / MPa")
plt.savefig("./stress/stress_loading_NGA.pdf", bbox_inches="tight", dpi=600)
#plt.show()
#%%

