from quickff.io import read_orca_hess_grad, make_yaff_ei
from quickff.log import log
from yaff import System
import numpy as np
import os
import h5py as h5
from quickff.tools import project_negative_freqs, get_ei_radii, average, charges_to_bcis

# import DFT data
numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_orca_hess_grad('dft_data/dut149_mol.input.hess', 'dft_data/dut149_mol.input.engrad')

# get atom types from Schmid group fragmentize protocol
atomtypes = np.genfromtxt('fragment_information/dut149_mol.input.mfpx',skip_header=2, dtype=str)[:,-2:]
ffatypes = []
for atom in atomtypes:
    ffatypes.append(atom[0]+"_"+atom[1])


with log.section('INIT', 4, timer='Initialization'):
    system = System(numbers, coords, rvecs, radii=None, masses=masses, ffatypes=ffatypes)
    system.detect_bonds()

with log.section('CHRG', 4, timer='Charges'):
    path='charges'
    with h5.File('dft_data/dut149_mol.h5', 'r') as f:
        charges = f[path][:]
        radii = None
        path_radii = os.path.join(os.path.dirname(path), 'radii')
        if 'radii' in f[path]:
            radii = average(f['%s/radii' %path][:], ffatypes, fmt='dict')
        else:
            radii = average(get_ei_radii(system.numbers), ffatypes, fmt='dict')

    bcis = charges_to_bcis(charges, ffatypes, system.bonds, verbose='True')
    make_yaff_ei('dut149_FF/dut149_FF_charges.out', None, bcis=bcis, radii=radii)