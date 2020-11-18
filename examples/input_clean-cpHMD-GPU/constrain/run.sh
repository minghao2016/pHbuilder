#!/bin/bash

# source specified gromacs version
source /usr/local/gromacs_test2/bin/GMXRC

gmx grompp -f MD.mdp -c 1cvo_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r 1cvo_NPT.pdb -maxwarn 1

gmx mdrun -v -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr -nb cpu
