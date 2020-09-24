#!/bin/python3

import os
from pythonlib import PDB

# Select MODEL 1 and remove header
protein = PDB()
protein.load("1cvo.pdb", 1)
protein.setTitle("CARDIOTOXIN V in water")
protein.write("1cvo_new.pdb", False)

model_ff    = "charmm36-mar2019"
model_water = "tip3p"

# Generate topology and protonate (make neutral) all GLU and ASP:
os.system("gmx pdb2gmx -f 1cvo_new.pdb -o 1cvo_new2.pdb -asp -glu -ignh -ff %s -water %s << EOF\n1\n1\n1\nEOF" % (model_ff, model_water))

# Add periodic box and put in center:
os.system("gmx editconf -f 1cvo_new2.pdb -o 1cvo_new3.pdb -c -d 1.0 -bt cubic")

# Add buffer waters to .pdb file:
os.system("gmx insert-molecules -f 1cvo_new3.pdb -o 1cvo_new4.pdb -ci buffer.pdb -nmol 3")

# Add the buffer water's topology to our .top file:
topList = []
with open("topol.top", "r") as file:
    for line in file.readlines():
        topList.append(line)

with open("topol.top", "w+") as file:
    for idx in range(0, len(topList)):
        file.write(topList[idx])

        if topList[idx] == "#endif\n" and topList[idx-1] == "#include \"posre.itp\"\n":
            file.write("\n; Include buffer topology\n")
            file.write("#include \"buffer.itp\"\n")

    file.write("BUF\t\t\t\t\t  %s\n" % (3))
topList.clear()

# Add normal waters:
os.system("gmx solvate -cp 1cvo_new4.pdb -o 1cvo_new5.pdb -p topol.top")

# Add ions:
os.system("gmx grompp -f ions.mdp -c 1cvo_new5.pdb -p topol.top -o ions.tpr")
os.system("gmx genion -s ions.tpr -o 1cvo_new6.pdb -p topol.top -pname NA -nname CL -neutral << EOF\nSOL\nEOF")

# Generate .ndx file














# pdbName     = "cardio"
# watermodel  = "spce"

# # Generate topology (9 = charmm36-mar2019.ff)
# os.system("echo 9 | gmx pdb2gmx -f %s.pdb -o %s_processed.pdb -asp -glu -ignh -water %s" % (pdbName, pdbName, watermodel))

# # gmx pdb2gmx -f cardio.pdb -o cardio_processed.pdb -asp -glu -ignh -water spce

# # Configure simulation box
# os.system("gmx editconf -f %s_processed.pdb -o %s_newbox.pdb -c -d 1.0 -bt cubic")

# 1. detect number of protonatable residues (ASP and GLU)

# 1. Fill box with water.
# 

# CHANGE .PDBtools FILE ############################################################

# Cardiotoxin:
#   By itself has a charge of 
#   62x residues
#   3x  H20 molecules of 3 atoms each (9 total) in buffer group
#   Bunch of water molecules
#   12x CL ion for neutral charge?

# CHANGE/GENERATE .TOP FILE
