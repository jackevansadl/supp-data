import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

plt.rcParams.update({'font.size': 6})

f, axarr = plt.subplots(3, sharex=True)
f.set_size_inches(1.8,2.6)
f.tight_layout()

data = np.loadtxt('porelocation.csv', delimiter=',')
axarr[0].fill_between(data[:,0], data[:,1], color="C1",label="cub.")
axarr[1].fill_between(data[:,0], data[:,2], color="C2",label="tet.")
axarr[2].fill_between(data[:,0], data[:,3], color="C4", label="oct.")
axarr[0].legend(loc=2,frameon=False)
axarr[1].legend(loc=2,frameon=False)
axarr[2].legend(loc=2,frameon=False)
axarr[2].set_xlabel("total loading / molecules$\,$UC$^{-1}$")
axarr[2].set_ylabel("loading in pore /  molecules$\,$UC$^{-1}$")
plt.savefig('./porecount.pdf', dpi=600, bbox='tight')

