
import sys
import numpy as np
import h5py as h5

from yaff import log
log.set_level(0)
from osmotic_mcmd.mcmd import MCMD
from molmod.units import kelvin, bar, kjmol, angstrom, femtosecond

system_file = "DUT-49_op.chk"
ff_file = "pars.txt"
adsorbate_file = "CH4.chk"

T = float(sys.argv[1]) * kelvin
P = float(sys.argv[2]) * bar
MD_trial_fraction = 0.00005
rcut = 12 * angstrom

from yaff.pes.eos import PREOS
eos = PREOS.from_name('methane')
fugacity = eos.calculate_fugacity(T,P)
mu = eos.calculate_mu(T,P)

mcmd = MCMD(system_file, adsorbate_file, ff_file, T, P, fugacity, MD_trial_fraction, rcut, write_h5s = False, barostat = True, vol_constraint = True)
mcmd.run_GCMC(40000000, 1000)



