#!/bin/python3

import phbuilder

# Set some pH-related data members in universe:
phbuilder.universe.add('ph_constantpH', True)
phbuilder.universe.add('ph_restrainpH', False)

phbuilder.universe.add('ph_GLU_dvdl', [26.238, -556.92, -106.76, 230.33, -155.89, -24.960])
phbuilder.universe.add('ph_ASP_dvdl', [44.936, -551.57, -109.62, 203.77, -127.44, -31.648])
phbuilder.universe.add('ph_BUF_dvdl', [672.41, -702.45, -63.10, 695.67, -1214.43, 537.14])

################################################################################

phbuilder.protein.process('../../proteins/ASP_tri.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019", "tip3p", d_terministring="34")

# syn-anti
phbuilder.topol.restrain_dihedrals('GLU', [' OE1', ' CD ', ' OE2', ' HE2'], 1,  0, 0, 10)
phbuilder.topol.restrain_dihedrals('ASP', [' OD1', ' CG ', ' OD2', ' HD2'], 1,  0, 0, 10)

# Ca-Cb
phbuilder.topol.restrain_dihedrals('GLU', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)
phbuilder.topol.restrain_dihedrals('ASP', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)

phbuilder.protein.add_box(d_boxMargin=1.0)
phbuilder.protein.add_buffer("../../proteins/buffer.pdb", "../../proteins/buffer.itp")
phbuilder.protein.add_water()
phbuilder.protein.add_ions()

phbuilder.md.energy_minimize()
phbuilder.md.energy_tcouple()
phbuilder.md.energy_pcouple()

phbuilder.md.gen_mdp('MD', nsteps=50000, nstxout=10000)
phbuilder.md.gen_constantpH(ph_pH=3.65, ph_lambdaM=5.0, ph_nstout=1, ph_barrierE=5.0)
phbuilder.write.run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")

phbuilder.write.jobscript('test', 48, 1, 32, 'lindahl')
phbuilder.universe.inspect()
