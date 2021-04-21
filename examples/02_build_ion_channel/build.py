#!/bin/python3

import phbuilder # Note: preparing this simulation might take a while.

# Set some pH-related data members in universe:
phbuilder.universe.add('ph_constantpH', True)
phbuilder.universe.add('ph_restrainpH', True)

phbuilder.universe.add('ph_GLU_dvdl', [25.65, -558.23, -201.97, 873.77, -1588.05, 1299.02, -436.91])
phbuilder.universe.add('ph_ASP_dvdl', [40.36, -549.89, -120.65, 283.59, -250.39, 38.67, -12.41])
phbuilder.universe.add('ph_BUF_dvdl', [673.29, -708.46, -160.24, 1395.34, -2797.63, 2017.18, -493.35])

################################################################################

phbuilder.protein.process("../../proteins/4hfi.pdb")

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019", "tip3p", d_terministring="11")

# syn-anti
phbuilder.topol.restrain_dihedrals('GLU', [' OE1', ' CD ', ' OE2', ' HE2'], 1,  0, 0, 10)
phbuilder.topol.restrain_dihedrals('ASP', [' OD1', ' CG ', ' OD2', ' HD2'], 1,  0, 0, 10)

# Ca-Cb
phbuilder.topol.restrain_dihedrals('GLU', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)
phbuilder.topol.restrain_dihedrals('ASP', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)

phbuilder.protein.add_box(d_boxMargin=1.5, d_boxType='triclinic')
phbuilder.protein.add_buffer("../../proteins/buffer.pdb", "../../proteins/buffer.itp", attempts=200000, d_bufMargin=1.5)
phbuilder.protein.add_water()
phbuilder.protein.add_ions()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()

phbuilder.md.gen_mdp('MD', nsteps=25000000, nstxout=10000)
phbuilder.md.gen_constantpH(ph_pH=4.0, ph_lambdaM=5.0, ph_nstout=200, ph_barrierE=7.5)
phbuilder.write.run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")

phbuilder.write.jobscript('test', 48, 1, 32, 'lindahl')
phbuilder.universe.inspect()
