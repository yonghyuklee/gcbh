import re
import os
import glob

from ase.io import read
from ase.io.trajectory import TrajectoryWriter
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator as SPC
from ase.constraints import FixAtoms
from mace.calculators import mace_mp
from ase.md.npt import NPT
from ase import units
from ase.md import MDLogger

def ramp_T(structure, calc, T0=1, T=300, ramp_steps=1000):
    """
        structure (ase structure) 
        calc (ase calculator)
        T0 - start temp
        T - end temp 
        ramp_steps - how many steps to reach ramp, split into 100 intervals
    """
    t = structure.copy()
    t.set_calculator(calc)
    
    ramp_step_temp = int(ramp_steps / 100)
    for T_temp in np.linspace(T0, T, 100):
        dyn = NPT(t, temperature_K=T_temp, externalstress=1., timestep=units.fs*2)    
        dyn.attach(MDLogger(dyn, t, 't_{}.log'.format(T_temp), header=False, stress=False,
           peratom=True, mode="a"), interval=3)   
        
        dyn.run(ramp_step_temp)
    
    return t


def ramp_P(structure, calc, T, P0=0, P=100, ramp_steps=1000):
    """
        structure (ase structure) 
        calc (ase calculator)
        P0 - start pressure
        P - end pressure 
        ramp_steps - how many steps to reach ramp, split into 100 intervals
    """
    t = structure.copy()
    t.set_calculator(calc)
    
    ramp_step_temp = int(ramp_steps / 100)
    for P_temp in np.linspace(P0, P, 100):
        dyn = NPT(t, temperature_K=T, externalstress=P_temp, timestep=units.fs*10)   
        dyn.attach(MDLogger(dyn, t, 'p_{}.log'.format(P_temp), header=False, stress=False,
           peratom=True, mode="a"), interval=3)        
        dyn.run(ramp_step_temp)
    return t


def md(structure, calc, T=300, P=0, steps=1000):
    """
        structure (ase structure) 
        calc (ase calculator)
        T - temp
        P - pressure 
        steps - number of steps
    """
    t = structure.copy()
    t.set_calculator(calc)
    # note the timestep here, it's 10 fs 
    dyn = NPT(t, temperature_K=T, externalstress=P, timestep=units.fs*10)    
    dyn.attach(MDLogger(dyn, t, 'md.log', header=False, stress=False,
           peratom=True, mode="a"), interval=3)
    dyn.run(steps)

    return t


def main():
    atoms = read("../../data/2f6l.xyz")
    vol=1
    atoms.set_cell((vol, vol, vol))
    atoms.center()
    calc = mace_mp(model="large", device="cuda")
    print("ramping Temp...")
    images = ramp_T(atoms, calc, T=500, ramp_steps=30000)
    print("ramping Pressure...")
    images = ramp_P(images, calc, T=500, P0=1, P=300, ramp_steps=30000)
    print("Running MD...")
    images = md(images, calc, T=300, P=300, timestep=1, steps=30000)
    images.write("opt.traj")


main()
