import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
fig = plt.gcf()
fig.tight_layout()
fig.set_size_inches(3.2,3)
plt.rcParams.update({'font.size': 6})

# data = pd.read_pickle("./dut49op_volume.pkl")
# #plt.plot(data['loading'], data['volume'], 'o', label='DUT-49op')
# plt.errorbar(data['loading'], data['volume'], yerr=data['volume_sd'], fmt='o',  label="op", color="C0")

data = pd.read_pickle("./dut49op_volume.pkl")
#plt.plot(data['loading'], data['volume'], 'o', label='DUT-49op')
plt.errorbar(data['loading'], (data['volume']-data['volume'].min())/data['volume'].min(), yerr=data['volume_sd']/data['volume'], fmt='o',  label="cp", color="C1")

# data = pd.read_pickle("./dut49cp_volume.pkl")
# plt.plot(data['loading'], data['volume'], 'o-', label='DUT-49cp')

plt.xlabel('loading / molecules$\,$UC$^{-1}$')
plt.ylabel('cell volume / Å$^3$')
plt.legend(loc=1)

plt.show()
#plt.savefig('./dut49volume_loading_cp.pdf', dpi=600, bbox_inches='tight')
