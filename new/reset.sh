#!/bin/bash

if [ -f "1cvo_MD.pdb" ]
then
	read -p "Warning: simulation has finished. Proceed? (y)" var
else
	rm -rf \_\_py* charmm*
	rm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg
	rm -f \#*\#
	rm -f buffer.pdb 1cvo_*.pdb
	rm -f run.sh reset.sh jobscript.sh universe
fi

if [ "${var}" = "y" ]
then
	rm -rf \_\_py* charmm*
	rm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg
	rm -f \#*\#
	rm -f buffer.pdb 1cvo_*.pdb
	rm -f run.sh reset.sh jobscript.sh universe
fi

