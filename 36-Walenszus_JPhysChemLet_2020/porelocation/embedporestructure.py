from ase.io import read, write
small = []
medium = []
large = []
adss = []
for ads in [10,20,30,40,50,60,70,80,90,100,120,130,140,150,160,170,180,190,200,250,350,400,450]:
    pore_count = []
    for unique in [1,2,3,4,5]:
        print(ads,unique)
        framework = read("porelocator/framework_"+str(ads)+"_"+str(unique)+".cif")
        adsorbate = read("porelocator/adsorbate_"+str(ads)+"_"+str(unique)+".cif")

        pores = read("porelocator/framework_"+str(ads)+"_"+str(unique)+"_porestructure.xyz")
        pores.set_cell(framework.cell)
        pores.set_pbc([1,1,1])

        combined = adsorbate+pores
        #write("combined.xyz", combined)
        # from ase.visualize import view
        # view(combined)

        CH2_sites = [atom.index for atom in adsorbate if atom.symbol == 'He']

        from ase.geometry import get_distances
        import numpy as np

        for i,k in zip(CH2_sites[0::2], CH2_sites[1::2]):
            com = combined[i,k].get_center_of_mass()
            distances = get_distances(com,pores.get_positions(),cell=pores.cell,pbc=[1,1,1])
            closest_pore = (np.argmin(distances[1]))
            pore_count.append(pores[closest_pore].number)

        del framework
        del pores
        del adsorbate

    small.append(float(pore_count.count(59))/5)
    medium.append(float(pore_count.count(60))/5)
    large.append(float(pore_count.count(61))/5)
    adss.append(ads)


print(small)
import matplotlib.pyplot as plt
plt.plot(adss,small, 'o-')
plt.plot(adss,medium, 'o-')
plt.plot(adss,large, 'o-')

import csv

with open('porelocation.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(adss, small,medium,large))
plt.show()