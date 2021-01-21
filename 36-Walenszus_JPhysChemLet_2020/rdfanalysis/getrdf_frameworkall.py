import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
fig = plt.gcf()
fig.tight_layout()
fig.set_size_inches(3.3,3)
plt.rcParams.update({'font.size': 6})

colors = plt.cm.Blues(np.linspace(0.3,1,24))
count = 0
for i in [10,20,30,40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 350, 400, 450]:
#for i in [10,100,200]:
    rdf1s = []
    rdf2s = []
    rdf3s = []
    rdf4s = []
    for j in [1,2,3,4,5]:
        with open("../rawdata/dut49/loading_"+str(i)+"_"+str(j)+"/tmp_"+str(j)+".rdf", "r") as f:
            lines = f.read().split('\n')

        rdf1 = []
        rdf2 = []
        rdf3 = []
        rdf4 = []
        for n in range(2000):
            nbutane1 = lines[4+(201*n):204+(201*n)]
            nbutane1 = [i.split() for i in nbutane1]
            nbutane1 = np.array(nbutane1, dtype=float)
            rad = (nbutane1[:,1])
            rdf1.append(nbutane1[:,12])
            rdf2.append(nbutane1[:,4])
            rdf3.append(nbutane1[:,6])
            rdf4.append(nbutane1[:,8])
        rdf1s.append(np.mean(rdf1, axis=0))
        rdf2s.append(np.mean(rdf2, axis=0))
        rdf3s.append(np.mean(rdf3, axis=0))
        rdf4s.append(np.mean(rdf4, axis=0))
        # steps = range(len(msd))

    mean_rdf1s = np.mean(rdf1s,axis=0)
    mean_rdf2s = np.mean(rdf2s,axis=0)
    mean_rdf3s = np.mean(rdf3s,axis=0)
    mean_rdf4s = np.mean(rdf4s,axis=0)


    # plt.plot(rad,np.mean(rdf3, axis=0), label="15-14")
    # plt.plot(rad,mean_rdf4s, label="CH$_2$-CH$_2$")
    plt.plot(rad,mean_rdf1s, color=colors[count], label=str(i))
    count=count+1

plt.xlabel('$r$ / Å')
plt.ylabel('g($r$)')

#plt.legend()
#plt.show()
plt.savefig('./rdf_nbutane_op.pdf', dpi=600, bbox='tight')

