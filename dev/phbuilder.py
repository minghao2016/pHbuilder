#!/usr/bin/python3
import os
from lib import sim

# PARAMETERS

gromPath    = "/home/anton/GIT/phbuilder/grom"
modelFF     = "charmm36-mar2019"
modelWater  = "tip3p"
pH          = 4.5           

# Copy files from grom to our working dir.
os.system("cp -r %s/* ." % gromPath)

################################################################################

sim = sim()
sim.processpdb("1cvo.pdb")

# PROTEIN
sim.protein_add_forcefield(modelFF, modelWater)
sim.protein_add_box()
sim.protein_add_buffer()
sim.protein_add_water()
sim.protein_add_ions()

# GENERATE INDEX FILE

group_PROTEIN     = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'SER', 'THR',
                     'VAL', 'PHE', 'TYR', 'TRP']
group_BUFFER      = ['BUF']
group_WATER       = ['SOL']
group_IONS        = [' NA', ' CL']
group_NON_PROTEIN = group_BUFFER + group_WATER + group_IONS
group_SYSTEM      = group_PROTEIN + group_NON_PROTEIN

sim.generate_index("SYSTEM"     , group_SYSTEM     )
sim.generate_index("PROTEIN"    , group_PROTEIN    )
sim.generate_index("BUFFER"     , group_BUFFER     )
sim.generate_index("WATER"      , group_WATER      )
sim.generate_index("IONS"       , group_IONS       )
sim.generate_index("NON_PROTEIN", group_NON_PROTEIN)

# GENERATE .MDP FILES

sim.generate_mdp("EM.mdp",  Type='EM',  dt=0.01,  nsteps=10000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate_mdp("NVT.mdp", Type='NVT', dt=0.002, nsteps=25000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate_mdp("NPT.mdp", Type='NPT', dt=0.002, nsteps=25000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate_mdp("MD.mdp",  Type='MD',  dt=0.002, nsteps=500000, output=1000, tgroups=[['SYSTEM', 0.5, 300]])

# GENERATE CONSTANT-PH .DAT FILE

sim.generate_phdata(pH)

# WRITE BASH SCRIPTS

sim.write_run("/usr/local/gromacs", "/usr/local/gromacs_dev")
sim.write_reset()
sim.write_jobscript("test", 36, 1)

# EQUILIBRATE

sim.energy_minimize()
sim.energy_tcouple()
sim.energy_pcouple()
