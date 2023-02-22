import sys
from quickff.io import read_orca_hess_grad, make_yaff_ei
from quickff.log import log
from yaff import System, ForceField
import numpy as np
import os
import h5py as h5
from quickff.tools import project_negative_freqs
from quickff.reference import SecondOrderTaylor, YaffForceField
from quickff.settings import Settings
from molmod import angstrom
from quickff.program import DeriveFF

log_out = 'dut49_FF/dut49_FF_quickFF.log'
log.set_level('highest')
log.line_length = 200
log.write_to_file(log_out)

# import DFT data
numbers, coords, energy, grad, hess, masses, rvecs, pbc = read_orca_hess_grad('dftdata/dut49_mol_again.input.hess', 'dftdata/dut49_mol_again.input.engrad')

# get atom types from Schmid group fragmentize protocol
atomtypes = np.genfromtxt('dftdata/dut49_mol_again.input.mfpx',skip_header=2, dtype=str)[:,-2:]
ffatypes = []
for atom in atomtypes:
    ffatypes.append(atom[0]+"_"+atom[1])

with log.section('INIT', 4, timer='Initialization'):
    system = System(numbers, coords, rvecs, radii=None, masses=masses, ffatypes=ffatypes)
    system.detect_bonds()
    hess = project_negative_freqs(hess, masses)

with log.section('REF', 4):
    energy=0.0
    ai = SecondOrderTaylor('ai', coords=coords, energy=energy, grad=grad, hess=hess)

fn_out = 'dut49_FF/dut49_FF_quickFF.txt'
chk_out = 'dut49_FF/dut49_ligand.chk'
settings = Settings(fn_yaff=fn_out, fn_sys=chk_out, plot_traj='Final', xyz_traj=False, consistent_cross_rvs=True, vdw="dut49_FF/dut49_FF_LJ.txt", vdw_rcut=20*angstrom, ei='dut49_FF/dut49_FF_charges.out', ei_rcut=50*angstrom, bond_term='bondmm3', bend_term="bendmm3", log_level='high', pert_traj_tol=1e-6)

with log.section('LJ_EI', 4):
    refs = []
    ff = ForceField.generate(system, settings.ei, rcut=settings.ei_rcut, alpha_scale=3.2, gcut_scale=1.0, tailcorrections=True)
    refs.append(YaffForceField('EI', ff))

    ff = ForceField.generate(system, settings.vdw, rcut=settings.vdw_rcut)
    refs.append(YaffForceField('vdW', ff))

#define quickff program
program = DeriveFF(system, ai, settings, ffrefs=refs)
#run program
program.run()
