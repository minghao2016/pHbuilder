#!/bin/bash

# Load constant-pH version of Gromacs
source /usr/local/gromacs_dev/bin/GMXRC

# Compile
gmx grompp -f MD.mdp -c 2khm_NPT.pdb -p topol.top -n index.ndx -o MD.tpr

# Run
gmx mdrun -v -s MD.tpr -o MD.trr -c 2khm_MD.pdb -g MD.log -e MD.edr

# Load normal version
source /usr/local/gromacs/bin/GMXRC
