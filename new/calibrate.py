#!/bin/python3

import phbuilder
import os, numpy, matplotlib.pyplot as plt, load

phbuilder.universe.add('d_constantpH', True)
phbuilder.universe.add('d_restrainpH', True)

phbuilder.universe.add('ph_GLU_dvdl', [26.238, -556.92, -106.76, 230.33, -155.89, -24.960])      # Noora new.
phbuilder.universe.add('ph_ASP_dvdl', [44.936, -551.57, -109.62, 203.77, -127.44, -31.648])      # Noora new.
# phbuilder.universe.add('ph_BUF_dvdl', [670.1, -674.4, 83.19, -150.21])                         # Noora original.
phbuilder.universe.add('ph_BUF_dvdl', [672.405, -702.445, -63.104, 695.686, -1214.433, 537.137]) # Calibrated using GLU_tri_capped with syn-anti and ca-cb.

phbuilder.protein.process('GLU_tri.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019", "tip3p", d_terministring="34")

# syn-anti
phbuilder.topol.restrain_dihedrals('GLU', [' OE1', ' CD ', ' OE2', ' HE2'], 1,  0, 0, 10)
phbuilder.topol.restrain_dihedrals('ASP', [' OD1', ' CG ', ' OD2', ' HD2'], 1,  0, 0, 10)

# Ca-Cb
phbuilder.topol.restrain_dihedrals('GLU', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)
phbuilder.topol.restrain_dihedrals('ASP', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)

phbuilder.protein.add_box(d_boxMargin=2.0)
phbuilder.protein.add_buffer("/home/anton/GIT/phbuilder/grom/buffer.pdb", "/home/anton/GIT/phbuilder/grom/buffer.itp", d_bufMargin=2.5)
phbuilder.protein.add_water()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()
phbuilder.write.run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")

# The part where we do the loop to get the mean and standard deviations:

dVdlInitList = [ii / 10.0 for ii in range(-1, 12)]
dVdlMeanList = []
dVdlStdList  = []

for init in dVdlInitList:
    phbuilder.md.energy_tcouple()
    phbuilder.md.energy_pcouple()
    phbuilder.md.gen_mdp('MD', nsteps=200000, nstxout=10000)
    phbuilder.md.gen_constantpH(pH=4.25, lambdaM=0.0, nstOut=1, barrierE=0.0, cal=True, lambdaInit=init)

    os.system("./run.sh")

    dVdlList = load.Col('lambda_2.dat', 3)
    dVdlMeanList.append(numpy.mean(dVdlList))
    dVdlStdList.append(numpy.std(dVdlList))

phbuilder.universe.add('d_dVdlInitList', dVdlInitList)
phbuilder.universe.add('d_dVdlMeanList', dVdlMeanList)
phbuilder.universe.add('d_dVdlStdList', dVdlStdList)

phbuilder.analyze.fitCalibration(5, compare=[670.1, -674.4, 83.19, -150.21])
