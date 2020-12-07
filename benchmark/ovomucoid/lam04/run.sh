#!/bin/bash

# This the "clean-cpHMD-GPU" branch gromacs version.
source /usr/local/gromacs_test2/bin/GMXRC

gmx grompp -f MD.mdp -c ovomucoid_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -maxwarn 1

# We run PME on the CPU and let gromacs decide for non-bonded and bonded (those are put on the GPU).
gmx mdrun -v -s MD.tpr -o MD.trr -c ovomucoid_MD.pdb -g MD.log -e MD.edr -pme cpu
