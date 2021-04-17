#!/bin/bash

# source specified gromacs version
source /usr/local/gromacs_test2/bin/GMXRC

gmx grompp -f MD.mdp -c protein_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r protein_NPT.pdb

gmx mdrun -v -s MD.tpr -o MD.trr -c protein_MD.pdb -g MD.log -e MD.edr -nb cpu
