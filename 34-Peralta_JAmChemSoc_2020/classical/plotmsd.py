import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("./20mol/1-butene/MSDOrderN/System_0/msd_self_1-butene_0.dat")

plt.plot(data[:,0],data[:,1], 'o-', label="1-butene")

data = np.genfromtxt("./20mol/cis-2-butene/MSDOrderN/System_0/msd_self_cis-2-butene_0.dat")

plt.plot(data[:,0],data[:,1], 'o-', label="cis-2-butene")

data = np.genfromtxt("./20mol/trans-2-butene/MSDOrderN/System_0/msd_self_trans-2-butene_0.dat")

plt.plot(data[:,0],data[:,1], 'o-', label="trans-2-butene")

plt.legend()
plt.xlabel("time / ps")
plt.ylabel("mean squared displacement / $\AA^2$")
plt.savefig("test.png")
