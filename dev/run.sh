#!/bin/bash

# NOORA ########################################################################

source /usr/local/gromacs_dev/bin/GMXRC

gmx grompp -f MD.mdp -c 1cvo_NPT.pdb -p topol.top -n index.ndx -o MD.tpr

gmx mdrun -v -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr

# PAUL #########################################################################

# source /usr/local/gromacs_paul/bin/GMXRC

# gmx grompp -f MD.mdp -c 1cvo_NPT.pdb -p topol.top -n index.ndx -o MD.tpr

# gmx mdrun -v -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr -cpHMD constant_ph_input.dat

################################################################################

# Load normal version
source /usr/local/gromacs/bin/GMXRC
