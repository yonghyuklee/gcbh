from ase.io import *
from ase.calculators.vasp import Vasp
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator as SPC
#from ase.constraints import FixAtoms

from ase import units
calc = Vasp(
     txt=      'stdout',
     prec=     'Normal',
     xc=       'PBE',
     pp=       'PBE',
     gamma=    True,
     reciprocal= True,
     ismear=   0,
     sigma=    0.1,
     nsw=      500,
     ibrion=   2,
     lwave=    False,
     lcharg=   False,
     lreal=    'Auto',
     lasph=    True,
     lorbit=   10,
     setups=   {'base': 'recommended'},
     ncore=    8,
     )

def vasp_relax(atoms, calc):
    t = atoms.copy()
    t.set_calculator(calc)
    t.get_potential_energy()
    return t

def main():
    atoms = read("./input.traj")
    kpts = [1,1,1]
    calc.set(kpts=kpts)
    images = vasp_relax(atoms, calc)
    e = images.get_potential_energy()
    f = images.get_forces()
    pos = images.get_positions()
    images.set_calculator(SPC(images, energy=e, forces=f))
    images.write("optimized.traj")

main()
