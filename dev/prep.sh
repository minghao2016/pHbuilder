#!/bin/bash

# GENERATE TOPOLOGY
# Note: this is protein-specific as we get prompted which starting charge we want
# for each GLU and ASP.
gmx pdb2gmx -f cardio.pdb -o cardio_processed.pdb -asp -glu -ignh -water tip3p << EOF
9
1
1
1
EOF

# CONFIGURE SIMULATION BOX / Fine
gmx editconf -f cardio_processed.pdb -o cardio_newbox.pdb -c -d 1.0 -bt cubic

# SOLVATE / Fine
gmx solvate -cp cardio_newbox.pdb -o cardio_solv.pdb -p topol.top

# IONS / Fine (use ions.mdp from tutorial)
gmx grompp -f ions.mdp -c cardio_solv.pdb -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o cardio_ions.pdb -p topol.top -pname NA -nname CL -neutral << EOF
SOL
EOF

# CHANGE THE FIRST (numAsp + numGlu) solvent molecules in cardio_ions.pdb to
# buffer particles.

# ENERGY MINIMIZATION
gmx grompp -f minim.mdp -c cardio_ions.pdb -p topol.top -o em.tpr
gmx mdrun -s em.tpr -c em.pdb
