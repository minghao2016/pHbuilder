#!/usr/bin/python3

import os
from lib import *

# PARAMETERS

gromPath    = "/home/anton/GIT/phbuilder/grom"
modelFF     = "charmm36-mar2019"
modelWater  = "tip3p"
pH          = 4.5           

sim = sim()
protein = PDB("1cvo.pdb")

# PATH AND FILE STUFF

# Copy files from grom to our working dir.
os.system("cp -r %s/* ." % gromPath)

# Backup previous builder.log file and index.ndx (if exists)
backupFile("builder.log")

################################################################################

# PROTEIN

pdbName = protein.fname(ext = 0)
# protein.resetResId()
protein.writepdb("%s_PR1.pdb" % (pdbName))

sim.loadpdb("1cvo.pdb")

sim.protein_add_forcefield(pdbName, modelFF, modelWater)
sim.protein_add_box(pdbName)
sim.protein_add_buffer(pdbName)
sim.protein_add_water(pdbName)      # Solvate
sim.protein_add_ions(pdbName)       # Genion

# GENERATE INDEX FILES

sim.loadpdb("%s_ION.pdb" % (pdbName))

group_PROTEIN     = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'SER', 'THR',
                     'VAL', 'PHE', 'TYR', 'TRP']
group_BUFFER      = ['BUF']
group_WATER       = ['SOL']
group_IONS        = [' NA', ' CL']
group_NON_PROTEIN = group_BUFFER + group_WATER + group_IONS
group_SYSTEM      = group_PROTEIN + group_NON_PROTEIN

sim.generate_index("index.ndx", "SYSTEM"     , group_SYSTEM      )
sim.generate_index("index.ndx", "PROTEIN"    , group_PROTEIN     )
sim.generate_index("index.ndx", "BUFFER"     , group_BUFFER      )
sim.generate_index("index.ndx", "WATER"      , group_WATER       )
sim.generate_index("index.ndx", "IONS"       , group_IONS        )
sim.generate_index("index.ndx", "NON_PROTEIN", group_NON_PROTEIN )

# GENERATE .MDP FILES

sim.generate_mdp("EM.mdp",  Type='EM',  dt=0.01,  nsteps=10000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate_mdp("NVT.mdp", Type='NVT', dt=0.002, nsteps=25000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate_mdp("NPT.mdp", Type='NPT', dt=0.002, nsteps=25000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate_mdp("MD.mdp",  Type='MD',  dt=0.002, nsteps=500000, output=1000, tgroups=[['SYSTEM', 0.5, 300]])

# GENERATE CONSTANT-PH .DAT FILE

sim.generate_phdata('%s_ION.pdb' % (pdbName), pH)

# WRITE BASH SCRIPTS

sim.write_run(pdbName, "/usr/local/gromacs", "/usr/local/gromacs_dev")
sim.write_reset(pdbName)
sim.write_jobscript(pdbName, "test", 36, 1)

# EQUILIBRATE

sim.energy_minimize(pdbName)
sim.energy_tcouple(pdbName)
sim.energy_pcouple(pdbName)
