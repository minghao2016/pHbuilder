#!/bin/bash

gmx mdrun -v -s MD.tpr -o MD.trr -c 1cvo_MD.pdb -g MD.log -e MD.edr
