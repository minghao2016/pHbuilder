#!/bin/python3

import phbuilder

phbuilder.universe.add('d_constantpH', True)
phbuilder.universe.add('d_restrainpH', True)

phbuilder.protein.process('3lr2.pdb')

phbuilder.write.reset()
phbuilder.topol.generate("charmm36-mar2019", "tip3p", d_terministring="11")
# phbuilder.topol.restrain_dihedrals('GLU', [' OE1', ' CD ', ' OE2', ' HE2'], 1,  0, 0, 10)
# phbuilder.topol.restrain_dihedrals('GLU', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)

phbuilder.protein.add_box(d_boxMargin=1.0)
phbuilder.protein.add_buffer("/home/anton/GIT/phbuilder/grom/buffer.pdb", "/home/anton/GIT/phbuilder/grom/buffer.itp")
# phbuilder.protein.add_water()
# phbuilder.protein.add_ions()

# phbuilder.md.energy_minimize()
# phbuilder.md.energy_tcouple()
# phbuilder.md.energy_pcouple()

# phbuilder.md.gen_mdp('MD', nsteps=10000000, nstxout=10000)
# phbuilder.md.gen_constantpH(7.0, 5.0, 100, barrierE=7.5)
# phbuilder.write.run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")

# phbuilder.write.jobscript('4hfi', 48, 1, 32, 'lindahl')

phbuilder.universe.inspect()