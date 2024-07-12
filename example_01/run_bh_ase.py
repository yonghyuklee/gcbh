##########################################################################################
# this is an example script to run the Grand Canonical Basin Hopping
# this self-contained script generates the executables to run basin hopping
# you still need to add the chemical potentials, input.traj, and the bh_options.json
##########################################################################################
# generalizations to the code such as general lammmps input file, etc. to come or whatever

import glob, json, argparse
import numpy as np
from ase.io import read
from gcbh2.scripts.gcbh2 import GrandCanonicalBasinHopping
# from pygcga2 import randomize, mirror_mutate, remove_H, add_H, rand_clustering
from pygcga2 import cluster_random_perturbation

# atom_elem_to_num = {"H": 1, "C": 6, "O": 8, "Al": 13, "Pt": 78}
atom_elem_to_num = {"Ti": 22, "Al": 13, "B": 5}
from ase.io import read
import numpy as np

def write_opt_vasp_file():
    # opt.py file
    with open("opt.py", "w") as f:
        f.write("from ase.io import *\n")
        f.write("from ase.calculators.vasp import Vasp\n")
        f.write("import numpy as np\n")
        f.write(
            "from ase.calculators.singlepoint import SinglePointCalculator as SPC\n"
        )
        f.write("#from ase.constraints import FixAtoms\n\n")


        f.write("from ase import units\n")
        f.write("calc = Vasp(\n")
        f.write("     txt=      'stdout',\n")
        f.write("     prec=     'Normal',\n")
        f.write("     xc=       'PBE',\n")
        f.write("     pp=       'PBE',\n")
        f.write("     gamma=    True,\n")
        f.write("     reciprocal= True,\n")
#        f.write("     encut=    500,\n")
#        f.write("     ediff=    1E-5,\n")
        f.write("     ismear=   0,\n")
        f.write("     sigma=    0.1,\n")
#        f.write("     nelm=     400,\n")
#        f.write("     nelmdl=   -10,\n")
        f.write("     nsw=      500,\n")
#        f.write("     ispin=    2,\n")
        f.write("     ibrion=   2,\n")
        f.write("     lwave=    False,\n")
        f.write("     lcharg=   False,\n")
        f.write("     lreal=    'Auto',\n")
        f.write("     lasph=    True,\n")
        f.write("     lorbit=   10,\n")
        f.write("     setups=   {'base': 'recommended'},\n")
        f.write("     ncore=    8,\n")
        f.write("     )\n\n")

        f.write("def vasp_relax(atoms, calc):\n")
        f.write('    t = atoms.copy()\n')
        f.write('    t.set_calculator(calc)\n')
        f.write('    t.get_potential_energy()\n')
        f.write('    return t\n\n')
        
        f.write("def main():\n")
        f.write('    atoms = read("./input.traj")\n')
        f.write('    kpts = [1,1,1]\n')
        f.write('    calc.set(kpts=kpts)\n')
        f.write('    images = vasp_relax(atoms, calc)\n')
        f.write('    e = images.get_potential_energy()\n')
        f.write('    f = images.get_forces()\n')
        f.write('    pos = images.get_positions()\n')
        f.write('    images.set_calculator(SPC(images, energy=e, forces=f))\n')
        f.write('    images.write("optimized.traj")\n\n')

        f.write("main()\n")

def write_opt_file(atom_order):
    # opt.py file
    with open("opt.py", "w") as f:
        f.write("import re\n")
        f.write("import os\n")
        f.write("import glob\n")
        f.write("\n")
        f.write("from ase.io import *\n")
        f.write("from ase.io.trajectory import TrajectoryWriter\n")
        f.write("import numpy as np\n")
        f.write(
            "from ase.calculators.singlepoint import SinglePointCalculator as SPC\n"
        )
        f.write("from ase.constraints import FixAtoms\n")
        #f.write("from pymatgen.io.lammps.data import LammpsData\n")
        f.write("from mace.calculators import mace_mp\n")
        f.write("from ase.md.velocitydistribution import MaxwellBoltzmannDistribution\n")
        f.write("from ase.md.verlet import VelocityVerlet\n")
        f.write("from ase import units\n")

        f.write("def md_mace_universal(surface, calc, T=500, timestep=1, steps=30):\n")
        f.write('    t = surface.copy()\n')
        f.write('    t.set_calculator(calc)\n')
        f.write('    dyn = MaxwellBoltzmannDistribution(t, temperature_K=T, force_temp=True)\n')
        #f.write('    dyn = VelocityVerlet(t, timestep * units.fs)\n')
        f.write('    dyn.run(steps)\n')
        f.write('    return t\n')

        
        f.write("def main():\n")
        f.write('    atoms = read("./input.traj")\n')
        f.write('    calc = mace_mp(model="large", device="cuda")\n')
        f.write("    images = md_mace_universal(atoms, calc, T=500, timestep=1, steps=30)\n")
        f.write('    images.write("opt.traj")\n')
        f.write("    e = images.get_potential_energy()\n")
        f.write("    f = images.get_forces()\n")
        f.write("    pos = images.get_positions()\n")
        f.write("    posz = pos[:, 2]\n")
        f.write("    ndx = np.where(posz < 5.5)[0]\n")
        f.write("    c = FixAtoms(ndx)\n")
        f.write("    images.set_constraint(c)\n")
        f.write("    images.set_calculator(SPC(images, energy=e, forces=f))\n")
        f.write('    images.write("optimized.traj")\n')
        f.write("main()\n")


def write_optimize_sh():
    with open("optimize.sh", "w") as f:
        f.write("pwd\n")
        f.write("cp ../../opt.py .\n")
        f.write("python opt.py\n")


def run_bh(options):
    filescopied = ["opt.py"]
    name = glob.glob(options["input_traj"])
    print(name)
    slab_clean = read(name[0])

    # this is the main part of the code
    bh_run = GrandCanonicalBasinHopping(
        temperature=options["temperature"],
        t_nve=options["t_nve"],
        atoms=slab_clean,
        bash_script="optimize.sh",
        files_to_copied=filescopied,
        restart=True,
        chemical_potential="chemical_potentials.dat",
    )

    dict_bonds = {
        "Ti-Ti": 3.2,
        "Ti-Al": 2.81,
        "Ti-B": 2.44,
        "Al-Al": 2.42,
        "Al-B": 2.05,
        "B-B": 1.68,
    }
    scalar_low = 0.5
    scalar_high = 5
    
    bond_range = {
                  tuple(bond.split('-')): [length * scalar_low, length * scalar_high]
                  for bond, length in dict_bonds.items()
                 }

    bh_run.add_modifier(
        cluster_random_perturbation, 
        name="cluster_random_perturbation", 
        dr_percent=2.0, 
        minimum_displacement=1.0, 
        bond_range=bond_range, 
        max_trial=50, 
        weight=2.0)
    n_steps = 4000
    bh_run.run(n_steps)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--options", type=str, default="./bh_options.json")
    args = parser.parse_args()
    with open(args.options) as f:
        options = json.load(f)

    #atom_order = options["atom_order"]
    write_opt_vasp_file()
    #write_lammps_input_file(atom_order=atom_order)
    write_optimize_sh()
    run_bh(options)


main()
