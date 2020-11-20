#!/usr/bin/python3

import os
from lib import sim

# Set up environment/external files.
gromPath = "/home/anton/GIT/phbuilder/grom"
os.system("cp {0}/buffer.itp {0}/buffer.pdb {0}/IONS.mdp ./".format(gromPath))

################################################################################

sim = sim()
sim.setconstantpH(True, restrain=False)
sim.processpdb("GlyGluGly.pdb", resetResId=True)

sim.write_reset()
                                                      # VI to have neutral terminii!
sim.protein_add_forcefield("charmm36-mar2019", "tip3p", neutralTermini=True)
sim.protein_add_box(boxSizeMargin=1.0)
sim.protein_add_buffer(minSep=1.5)
sim.protein_add_water()
sim.protein_add_ions()

sim.generate_index()

sim.generate_mdp('EM')
sim.generate_mdp('NVT')
sim.generate_mdp('NPT') # 10'000 steps enough for obtaining good dV/dl average. (?)
sim.generate_mdp('MD', nsteps=10000, nstxout=0, nstvout=0)

sim.generate_phdata_legacy(4.5, lambdaM=0, nstOut=1, barrierE=5.0)

sim.write_run("/usr/local/gromacs_test", mode='cpu')
sim.write_reset()

sim.energy_minimize()
sim.energy_tcouple()
sim.energy_pcouple(skip=True)
