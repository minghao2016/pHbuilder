#!/bin/bash

# NOORA ########################################################################

# source /usr/local/gromacs_dev/bin/GMXRC

# gmx grompp -f MD.mdp -c 1cvo_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r 1cvo_NPT.pdb

# gmx mdrun -v -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr

# PAUL #########################################################################

source /usr/local/gromacs_paul/bin/GMXRC

# gmx grompp -f MD.mdp -c 1cvo_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r 1cvo_NPT.pdb

gmx mdrun -v -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr -cpHMD constant_ph_input.dat -nb cpu

################################################################################

# CONTINUE
# gmx convert-tpr -s MD.tpr -o MD.tpr -extend 2000
# gmx mdrun -v -cpi state.cpt -append -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr

# Load normal version
source /usr/local/gromacs/bin/GMXRC