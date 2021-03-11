
import numpy as np
np.random.seed(52) # Ensure each process starts from the same configuration
import sys, os
import h5py as h5
from molmod.units import kelvin, bar, kjmol, angstrom, femtosecond

from yaff import System, log
from yaff.sampling.io import HDF5Writer, XYZWriter
from yaff.sampling.nvt import *
from yaff.sampling.npt import *
from yaff.sampling.verlet import *
from yaff.pes.ff import ForceField
from yaff.external.lammps_generator import *
from yaff.external.lammpsio import *
from yaff.external.liblammps import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank==0:
    log.set_level(log.medium)
else:
    log.set_level(log.silent)

T = float(sys.argv[1])
P = float(sys.argv[2])
vol_constraint = sys.argv[3] == 'True'

def run():

    # Input files
    fn_host = 'DUT-49_op.chk'
    host = System.from_file(fn_host).supercell(1,1,1)
    fn_pars = 'pars.txt'

    if rank == 0:
        fh5 = h5.File('trajectory.h5','w')
        hdf5writer = HDF5Writer(fh5, step=100)
        vsl = VerletScreenLog(step=100)

    ff = ForceField.generate(host, fn_pars, rcut=12*angstrom, alpha_scale=3.2, gcut_scale=1.0, tailcorrections=True)

    ff_lammps = swap_noncovalent_lammps(ff, fn_system='system.dat',
                                    fn_table='table.dat',
                                    nrows=5000,
                                    kspace='pppm',
                                    kspace_accuracy=1e-7,
                                    move_central_cell=False,
                                    fn_log='none',
                                    overwrite_table=False, comm=comm)

    mtk = MTKBarostat(ff_lammps, temp=T, press=P, timecon=1000*femtosecond, vol_constraint = vol_constraint, anisotropic = True)
    nhc = NHCThermostat(temp=T, timecon=100*femtosecond, chainlength=3)
    tbc = TBCombination(nhc, mtk)

    if rank == 0:
        verlet = VerletIntegrator(ff_lammps, 0.5*femtosecond, hooks=[tbc, hdf5writer, vsl])
    else:
        verlet = VerletIntegrator(ff_lammps, 0.5*femtosecond, hooks=[tbc])
    verlet.run(2000000)

    if rank == 0:
        fh5.close()


run()

