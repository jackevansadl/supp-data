import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib as mpl
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42
fig = plt.gcf()
fig.tight_layout()
fig.set_size_inches(3.3,3)
plt.rcParams.update({'font.size': 6})

loading = []
diff = []
std = []
for i in [10,20,30,40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 350, 400, 450]:
    msds = []
    slopes = []
    for j in [1,2,3,4,5]:
        data = np.genfromtxt("../rawdata/dut49/loading_"+str(i)+"_"+str(j)+"/msd_avetime_"+str(j)+".out")
        msd_total = data[int(len(data)/3):,4]
        steps_total = data[int(len(data)/3):,0]
        slope, intercept, r_value, p_value, std_err = linregress(steps_total, msd_total)
        slopes.append((slope/6)*1e-5)
        msds.append(msd_total)
    diff.append(np.mean(slopes))
    std.append(np.std(slopes))
    loading.append(i)


print(diff)
print(std)
plt.errorbar(loading, diff, yerr=std, fmt='o',  label="op", color="k",markersize=4,elinewidth=1)
dataset = pd.DataFrame({'loading': loading, 'diffusion':diff, 'diffusion_variance':std})
dataset.to_csv("./diffusion/diffusion_op.csv")

loading = []
diff = []
std = []
for i in [10,20,30,40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]:
    msds = []
    slopes = []
    for j in [1,2,3,4,5]:
        data = np.genfromtxt("../rawdata/dut49_cp/loading_"+str(i)+"_"+str(j)+"/msd_avetime_"+str(j)+".out")
        msd_total = data[int(len(data)/3):,4]
        steps_total = data[int(len(data)/3):,0]
        slope, intercept, r_value, p_value, std_err = linregress(steps_total, msd_total)
        slopes.append((slope/6)*1e-5)
        msds.append(msd_total)
    diff.append(np.mean(slopes))
    std.append(np.std(slopes))
    loading.append(i)


print(diff)
print(std)
plt.errorbar(loading, diff, yerr=std, fmt='s',  label="cp", color="grey",markersize=4,elinewidth=1)
dataset = pd.DataFrame({'loading': loading, 'diffusion':diff, 'diffusion_variance':std})
dataset.to_csv("./diffusion/diffusion_cp.csv")

plt.legend(loc=4,frameon=False)
plt.ylabel('$D$ / m$^2\,$s$^{-1}$')
plt.xlabel("loading / molecules$\,$UC$^{-1}$")
plt.yscale('log')
plt.xlim([0, 250])
plt.show()
