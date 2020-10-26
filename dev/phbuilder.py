#!/usr/bin/python3

import os, sys
from lib import sim

# Set up environment/external files.
gromPath = "/home/anton/GIT/phbuilder/grom"
os.system("cp {0}/buffer.itp {0}/buffer.pdb {0}/IONS.mdp ./".format(gromPath))

################################################################################

sim = sim()
sim.setconstantpH(True)
sim.processpdb("1cvo.pdb")

# PROTEIN
sim.protein_add_forcefield("charmm36-mar2019", "tip3p")
sim.protein_add_box(boxSizeMargin=1.0)
sim.protein_add_buffer(minSep=1.5)
sim.protein_add_water()
sim.protein_add_ions()

# GROUPS
group_PROTEIN     = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'SER', 'THR',
                     'VAL', 'PHE', 'TYR', 'TRP']
group_BUFFER      = ['BUF']
group_WATER       = ['SOL']
group_IONS        = [' NA', ' CL']
group_NON_PROTEIN = group_BUFFER + group_WATER + group_IONS
group_SYSTEM      = group_PROTEIN + group_NON_PROTEIN

# GENERATE INDEX FILE
sim.generate_index("SYSTEM"     , group_SYSTEM     )
sim.generate_index("PROTEIN"    , group_PROTEIN    )
sim.generate_index("BUFFER"     , group_BUFFER     )
sim.generate_index("WATER"      , group_WATER      )
sim.generate_index("IONS"       , group_IONS       )
sim.generate_index("NON_PROTEIN", group_NON_PROTEIN)

# GENERATE .MDP FILES
sim.generate_mdp('EM')
sim.generate_mdp('NVT')
sim.generate_mdp('NPT')
sim.generate_mdp('MD',  nsteps=500000, output=1)

# GENERATE CONSTANT-PH .DAT FILE
sim.generate_phdata(pH=4.5, lambdaM=3, nstOut=1, barrierE=7.5)

# WRITE BASH SCRIPTS
sim.write_run("/usr/local/gromacs_dev")
sim.write_jobscript(jobName="test", time=1, nodes=1, ntasks=32, queue="tcb")
sim.write_reset()

# EQUILIBRATE
sim.energy_minimize()
sim.energy_tcouple()
sim.energy_pcouple()
