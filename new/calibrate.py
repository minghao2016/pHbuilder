#!/bin/python3

import phbuilder
import os, numpy, matplotlib.pyplot as plt, load

phbuilder.universe.add('d_constantpH', True)
phbuilder.universe.add('d_restrainpH', False)

phbuilder.protein.process('GLU_tri.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019", "tip3p", d_terministring="34")
phbuilder.topol.restrain_dihedrals('GLU', [' OE1', ' CD ', ' OE2', ' HE2'], 1,  0, 0, 10)
phbuilder.topol.restrain_dihedrals('GLU', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)

phbuilder.protein.add_box(d_boxMargin=1.0)
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
    phbuilder.md.gen_mdp('MD', nsteps=100000, nstxout=0)
    phbuilder.md.gen_constantpH(pH=4.25, lambdaM=0.0, nstOut=1, barrierE=7.5, cal=True, lambdaInit=init)

    os.system("./run.sh")

    dVdlList = load.Col('lambda_1.dat', 3)
    dVdlMeanList.append(numpy.mean(dVdlList))
    dVdlStdList.append(numpy.std(dVdlList))

phbuilder.universe.add('d_dVdlInitList', dVdlInitList)
phbuilder.universe.add('d_dVdlMeanList', dVdlMeanList)
phbuilder.universe.add('d_dVdlStdList', dVdlStdList)

phbuilder.analyze.fitCalibration(5, compare=[26.238, -556.92, -106.76, 230.33, -155.89, -24.96])
