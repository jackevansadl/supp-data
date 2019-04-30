import sys
import json
import yaml
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pprint import pprint
from sklearn.preprocessing import RobustScaler
import pymatgen as pmg
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.core.structure import Structure
import scipy.stats as stats
from pymatgen.util.coord_utils import pbc_shortest_vectors, get_angle
from pymatgen.io.zeopp import ZeoCssr, ZeoVoronoiXYZ, get_voronoi_nodes, \
    get_high_accuracy_voronoi_nodes, get_void_volume_surfarea, \
    get_free_sphere_params
from sklearn.ensemble import AdaBoostRegressor, GradientBoostingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error
import matplotlib as mpl
import itertools
from matminer.featurizers.structure import DensityFeatures, \
    RadialDistributionFunction, \
    RadialDistributionFunctionPeaks, PartialRadialDistributionFunction

import re

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["ps.fonttype"] = 42
SMALL_SIZE = 8
mpl.rc('font', size=SMALL_SIZE)

f, axarr = plt.subplots(3,3, sharey=True, sharex=True)
structure_bare = Structure.from_file("DUT67_desolvated.cif")

structure_dosed = Structure.from_file("DUT67_loading1.cif")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Kr":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Kr":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[0,0].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[0,0].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")


Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Xe":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Xe":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[0,1].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[0,1].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Rn":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Rn":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[0,2].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[0,2].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")
axarr[0,2].legend(frameon=False,loc=1)

#####################################################################################
structure_dosed = Structure.from_file("DUT67_loading2.cif")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Kr":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Kr":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[1,0].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[1,0].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Xe":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Xe":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[1,1].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[1,1].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Rn":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Rn":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[1,2].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[1,2].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")

#####################################################################################
structure_dosed = Structure.from_file("DUT67_loading3.cif")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Kr":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Kr":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[2,0].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[2,0].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Xe":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Xe":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[2,1].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[2,1].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")

Kr_index_bare = []
for index, site in enumerate(structure_bare):
    if site.species_string == "Rn":
        Kr_index_bare.append(index)

distances_bare = []
amount_bare = []
for i in Kr_index_bare:
     neighbors_bare = structure_bare.get_neighbors(structure_bare[i], 20, include_index=True)
     
     framework_list_bare = []
     for element in neighbors_bare:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_bare.append(element)
             

     for element in framework_list_bare:
         distances_bare.append(element[1])
         occupancy = str(structure_bare[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_bare.append(float(occupancy))

Kr_index_dosed = []
for index, site in enumerate(structure_dosed):
    if site.species_string == "Rn":
        Kr_index_dosed.append(index)

distances_dosed = []
amount_dosed = []
for i in Kr_index_dosed:
     neighbors_dosed = structure_dosed.get_neighbors(structure_dosed[i], 20, include_index=True)
             
     framework_list_dosed = []
     for element in neighbors_dosed:
         if element[0].species_string[:1] not in ["Kr", "Xe", "Rn"]:
             framework_list_dosed.append(element)


     for element in framework_list_dosed:
         distances_dosed.append(element[1])
         occupancy = str(structure_dosed[element[2]]._species)
         occupancy = re.sub("[^0123456789\.]","",occupancy)
         amount_dosed.append(float(occupancy))


dist_hist_bare, dist_bins_bare = np.histogram(distances_bare, bins=np.arange(0, 15 + 1, 1), weights=amount_bare, density=False)
dist_hist_dosed, dist_bins_dosed = np.histogram(distances_dosed, bins=np.arange(0, 15 + 1, 1), weights=amount_dosed, density=False)

axarr[2,2].plot(dist_bins_bare[:-1], dist_hist_bare, color="C0", label="bare")
axarr[2,2].plot(dist_bins_dosed[:-1], dist_hist_dosed, color="C1", label="dosed")


axarr[2,1].set_xlabel("radius / $\AA$")
axarr[1,0].set_ylabel("frequency")
# plt.show()
plt.gcf().set_size_inches(7,7)
plt.savefig("rff_dut67.pdf", format="pdf", bbox_inches="tight", dpi=600, transparent=True)
