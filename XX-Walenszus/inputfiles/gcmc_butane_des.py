
import numpy as np
import h5py as h5
import sys

from molmod.units import kelvin, bar, angstrom, kjmol, liter, pascal, amu
from molmod.constants import avogadro

from yaff.pes.eos import PREOS
from yaff.sampling.io import MCHDF5Writer
from yaff.sampling.mc import GCMC
from yaff.sampling.mcutils import MCScreenLog
from yaff.system import System
from yaff import log

T = 298

fn_guest = str(sys.argv[1])
fn_host = str(sys.argv[2])
fn_pars = [str(sys.argv[3]), str(sys.argv[4])]

P = float(sys.argv[5])*bar

log.set_level(log.medium)

def simulate():
    
    host = System.from_file(fn_host).supercell(1,1,1)

    mc_moves =  {'insertion':1.0, 'deletion':1.0,
                 'translation':1.0, 'rotation':1.0}

    # Construct equation of state to link pressure, fugacity and chemical potential for butane
    eos = PREOS(425.125*kelvin, 3796000.0*pascal, 0.201, mass=58.12*amu)
  
    fugacity = eos.calculate_fugacity(T,P)
    mu = eos.calculate_mu(T,P)
    # Screen logger
    screenlog = MCScreenLog(step=10000)
    # HDF5 logger
    fh5 = h5.File('trajectory_'+str(sys.argv[5])+'.h5','w')
    hdf5writer = MCHDF5Writer(fh5, step=10000)
    # Setup the GCMC calculation, this generates a lot of output so when
    # force fields are generated, so we silence the logging for a while.
    log.set_level(log.silent)
    gcmc = GCMC.from_files(fn_guest, fn_pars, host=host,
        rcut=12.0*angstrom, tr=None, tailcorrections=True, hooks=[screenlog, hdf5writer],
        reci_ei='ewald_interaction', nguests=600)
    log.set_level(log.medium)
    # Set the external conditions
    gcmc.set_external_conditions(T, fugacity)
    # Run MC simulation
    init_h5 = h5.File('high_pressure_state.h5','r')
    init_system = System.from_hdf5(init_h5['/snapshots/000010000000'])
    gcmc.run(int(1e7), mc_moves=mc_moves, initial=init_system)
    fh5.close()

if __name__=='__main__':
    simulate()
