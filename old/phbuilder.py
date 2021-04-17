#!/bin/python3

import os
from lib import sim

# Set up environment/external files.
gromPath = "/home/anton/GIT/phbuilder/grom"
os.system("cp {0}/buffer.itp {0}/buffer.pdb {0}/IONS.mdp ./".format(gromPath))

################################################################################

sim = sim()
sim.setconstantpH(False, restrain=False)
sim.processpdb("1cvo.pdb")

sim.protein_add_forcefield("charmm36-mar2019", "tip3p")
sim.protein_add_box(boxSizeMargin=1.5)
# sim.protein_add_buffer(minSep=2.0)
sim.protein_add_water()
sim.protein_add_ions()

sim.generate_index()

sim.generate_mdp('EM')
sim.generate_mdp('NVT')
sim.generate_mdp('NPT')
sim.generate_mdp('MD', nsteps=100000, nstxout=0)

sim.generate_phdata(pH=4.5, lambdaM=5.0, nstOut=5000, barrierE=5.0)
# sim.generate_phdata_legacy(pH=4.25, lambdaM=5.0, nstOut=1, barrierE=5.0)

# sim.write_jobscript('test', 1, 1, 32, 'lindahl')
sim.write_run("/usr/local/gromacs_test2", mode='gpu')
sim.write_reset()

sim.energy_minimize()
sim.energy_tcouple()
sim.energy_pcouple(skip=True)