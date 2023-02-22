
from yaff import System, ForceField, log
from molmod import angstrom
import h5py as h5
import sys
from yaff import ff2lammps, write_lammps_table
import os

log.set_level(log.medium)

system_fn = str(sys.argv[1])
ff_fn = [str(sys.argv[2]), str(sys.argv[3])]

name = system_fn.split("/")[-1]
name = name.replace(".chk","")

dn = './'

system = System.from_file(system_fn)
ff = ForceField.generate(system, ff_fn, rcut=12*angstrom, alpha_scale=3.2, gcut_scale=1.0, tailcorrections=True)
write_lammps_table(ff, fn='lammps.table', rmin=0.50*angstrom, nrows=5000, unit_style='real')
ff2lammps(system, ff_fn, dn, rcut=12.0*angstrom, tailcorrections=False,
            tabulated=True, unit_style='real')

with open(os.path.join(dn,'lammps.in'),'r') as f:
    lines = f.readlines()
with open(os.path.join(dn,'lammps.in'),'w') as f:
    for line in lines[:-5]:
        f.write(line)
    f.write("thermo 1000\ntimestep 0.5 # in time units\n")
    f.write("velocity all create 300.0 666 # initial temperature in Kelvin and random seed\n")
    f.write("fix 1 all npt temp 300.0 300.0 100.0 tri 1.0 1.0 1000.0 tchain 3 mtk yes nreset 1000\ndump      traj_dcd all dcd 1000 traj.dcd\n")
    f.write("fix_modify 1 energy yes # Add thermo/barostat contributions to energy\n")
    f.write("run 4000000\n")


