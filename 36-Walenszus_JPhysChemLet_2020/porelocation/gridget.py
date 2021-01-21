import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

colors = plt.cm.viridis(np.linspace(0.3,1,24))
for index,ads in enumerate([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]):
    for unique in [1,2,3,4,5]:
        data = np.genfromtxt("porelocator_cp/framework_"+str(ads)+"_"+str(unique)+".vpsdpts",skip_header=1)
      
        from sklearn.neighbors import KernelDensity
        import statsmodels.api as sm
        X = data[:,4].reshape(-1, 1)
        kde = sm.nonparametric.KDEUnivariate(X)
        kde.fit(adjust=2.0)
        plt.plot(kde.support, kde.density, color=colors[index])


        from scipy.signal import argrelextrema
        mi, ma = argrelextrema(kde.density, np.less)[0], argrelextrema(kde.density, np.greater)[0]

        print(kde.support[mi])
        diffs = np.diff(kde.support[mi])

        pores = []
        points = []
        print(len(mi))
        for point in range(len(X)):
            points.append([data[point,1],data[point,2],data[point,3]])
            if len(mi) == 2:
                if data[point,4] <= kde.support[mi[0]]:
                    pore = 'Ce'
                if kde.support[mi[0]] < data[point,4] <= kde.support[mi[1]]:
                    pore = 'Pr'
                if kde.support[mi[1]] < data[point,4]:
                    pore = "Nd"
            if len(mi) == 3:
                pore = None
                if data[point,4] <= kde.support[mi[0]]:
                    pore = 'Ce'
                if kde.support[mi[0]] < data[point,4] <= kde.support[mi[1]]:
                    pore = 'Pr'
                if kde.support[mi[1]] < data[point,4] <= kde.support[mi[2]]:
                    pore = "Nd"
                if kde.support[mi[2]] < data[point,4]:
                    pore = "Pm"
            if len(mi) == 4:
                pore = None
                if data[point,4] <= kde.support[mi[0]]:
                    pore = 'Ce'
                if kde.support[mi[0]] < data[point,4] <= kde.support[mi[1]]:
                    pore = 'Pr'
                if kde.support[mi[1]] < data[point,4] <= kde.support[mi[2]]:
                    pore = "Nd"
                if kde.support[mi[2]] < data[point,4] <= kde.support[mi[3]]:
                    pore = "Pm"
                if kde.support[mi[3]] < data[point,4]:
                    pore = "Sm"
            if len(mi) == 5:
                pore = None
                if data[point,4] <= kde.support[mi[0]]:
                    pore = 'He'
                if kde.support[mi[0]] < data[point,4] <= kde.support[mi[1]]:
                    pore = 'Ne'
                if kde.support[mi[1]] < data[point,4] <= kde.support[mi[2]]:
                    pore = "Ar"
                if kde.support[mi[2]] < data[point,4] <= kde.support[mi[3]]:
                    pore = "Kr"
                if kde.support[mi[3]] < data[point,4] <= kde.support[mi[4]]:
                    pore = "Xe"
                if kde.support[mi[4]] < data[point,4]:
                    pore = "Rn"
            if len(mi) == 6:
                pore = None
                if data[point,4] <= kde.support[mi[0]]:
                    pore = 'He'
                if kde.support[mi[0]] < data[point,4] <= kde.support[mi[1]]:
                    pore = 'Ne'
                if kde.support[mi[1]] < data[point,4] <= kde.support[mi[2]]:
                    pore = "Ar"
                if kde.support[mi[2]] < data[point,4] <= kde.support[mi[3]]:
                    pore = "Kr"
                if kde.support[mi[3]] < data[point,4] <= kde.support[mi[4]]:
                    pore = "Xe"
                if kde.support[mi[4]] < data[point,4] <= kde.support[mi[5]]:
                    pore = "Rn"
                if kde.support[mi[5]] < data[point,4]:
                    pore = "La"
            if len(mi) == 7:
                pore = None
                if data[point,4] <= kde.support[mi[0]]:
                    pore = 'He'
                if kde.support[mi[0]] < data[point,4] <= kde.support[mi[1]]:
                    pore = 'Ne'
                if kde.support[mi[1]] < data[point,4] <= kde.support[mi[2]]:
                    pore = "Ar"
                if kde.support[mi[2]] < data[point,4] <= kde.support[mi[3]]:
                    pore = "Kr"
                if kde.support[mi[3]] < data[point,4] <= kde.support[mi[4]]:
                    pore = "Xe"
                if kde.support[mi[4]] < data[point,4] <= kde.support[mi[5]]:
                    pore = "Rn"
                if kde.support[mi[5]] < data[point,4] <= kde.support[mi[6]]:
                    pore = "La"
                if kde.support[mi[6]] < data[point,4]:
                    pore = "Ce"
            pores.append(pore)

        from ase import Atoms
        test = Atoms(pores, positions=points)
        from ase.io import read, write
        write("porelocator/framework_"+str(ads)+"_"+str(unique)+"_porestructure.xyz", test)

plt.show()

