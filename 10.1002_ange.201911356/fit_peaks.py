import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
from scipy.optimize import leastsq

cmap = mpl.cm.get_cmap('tab20').colors
print(cmap[0])


def _2Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2))+(amp2*wid2**2/((x-cen2)**2+wid2**2))


#import data
data = np.genfromtxt("./data/PPP_MMM.csv", delimiter=";")

idx_hi = np.abs(data[:,0] - 1150).argmin()
idx_lo = np.abs(data[:,0] - 1350).argmin()


#crop data
wavenumbers = data[idx_lo:idx_hi,0]
intensities = data[idx_lo:idx_hi,1:]

# peak1 = 1225
# peak2 = 1280

def peakfit(xdata, ydata, peak1, peak2):
    
    #initialize values
    amp1 = 1
    wid1 = 1
    cen1 = peak1 
    amp2 = 1
    wid2 = 1
    cen2 = peak2

    xData, yData = xdata, ydata

    #normalize
    yData = yData / max(yData)

    #fit function
    popt_2Lorentzian, pcov_2Lorentzian = scipy.optimize.curve_fit(_2Lorentzian, xData, yData, p0=[amp1, cen1, wid1, amp2,cen2,wid2])


    return(popt_2Lorentzian)

names = ["PPP-MM-Ref-1", "PPP-MM-Ref-2", "PPP-MM-1", "PPP-MM-2", "PPP-MM-3", "PPP-MM-4", "PPP-MM-5","PPP-MM-3_6h", "PPP-MM-3_4h", "PPP-MM-3_2h", "PPP-PBM-11", "PPP-PBM-12", "PPP-PBM-14"]
amplitude1 = []
amplitude2 = []
peaks1 = []
peaks2 = []
peak_differences = []
peak_ratio = []
for i in range(len(names)):
    plt.plot(wavenumbers, intensities[:,i], 'o', color="k")
    results = peakfit(wavenumbers,intensities[:,i], 1225, 1280)

    amplitude1.append(results[0])
    amplitude2.append(results[3])
    peaks1.append(results[1])
    peaks2.append(results[4])

    peak_differences.append(results[4]-results[1])
    peak_ratio.append(results[3]/results[0])

    #names.append("PPP-MMM-"+str(i+1))

    plt.plot(wavenumbers, _2Lorentzian(wavenumbers, *results)*max(intensities[:,i]), '-', color=(cmap[i]))
    plt.xlabel("wavenumber / cm$^{-1}$")
    plt.ylabel("intensity")
    plt.title(names[i], loc="left")
    plt.savefig(names[i]+".pdf", dpi=600, transparent=True, bbox_inches='tight')
    plt.show()
    plt.close()

import pandas as pd

d = {'name':names, 'peak1':peaks1, 'amplitude1':amplitude1, 'peak2':peaks2, 'amplitude2':amplitude2,  'peak_difference':peak_differences, "peak_ratio":peak_ratio}
df = pd.DataFrame(d)

df.to_csv('PPP_MM_peakfits.csv', index=False)
