#!/bin/python3

import os
from lib import PDB

protein = PDB()

# PARAMS #######################################################################

pdbName    = "1cvo"
pdbDescrip = "CARDIOTOXIN V in water"

gromPath   = "/home/anton/GIT/phbuilder/grom"   # relative path to grom dir
modelFF    = "charmm36-mar2019"
modelWater = "tip3p"

startProtonated = True                  # ASP and GLU will be neutral

protein.loadpdb("%s.pdb" % (pdbName))

# PATH STUFF ###################################################################

# Create symbolic link to force field dir for GROMACS (cause it's large).
if os.path.islink("%s.ff" % (modelFF)):
    os.remove("%s.ff" % modelFF)

os.symlink("%s/%s.ff" % (gromPath, modelFF), "%s.ff" % (modelFF))

# Copy files from grom to our working dir.
os.system("cp %s/* ." % (gromPath))

# GMX PDB2GMX ##################################################################

protein.setTitle(pdbDescrip)
protein.writepdb("%s_PR1.pdb" % (pdbName))

countASP  = protein.countRes("ASP")
countGLU  = protein.countRes("GLU")
countACID = countASP + countGLU

# Print how many acidic residues were found
print("Detected %s acidic residues (%s ASP and %s GLU)..." % (countACID, countASP, countGLU))

# Create EOF string required for pdb2gmx
xstr = "<< EOF"
for idx in range(0, countACID):
    xstr += "\n%s" % (int(startProtonated))
xstr += "\nEOF"

# Generate topology and protonate (make neutral) all GLU and ASP:
os.system("gmx pdb2gmx -f %s_PR1.pdb -o %s_PR2.pdb -asp -glu -ignh -ff %s -water %s %s" % (
    pdbName, pdbName, modelFF, modelWater, xstr))

# GMX EDITCONF #################################################################

os.system("gmx editconf -f %s_PR2.pdb -o %s_BOX.pdb -c -d 1.0 -bt cubic" % (pdbName, pdbName))

# GMX INSERT-MOLECULES (ADD BUFFER) ############################################

os.system("gmx insert-molecules -f %s_BOX.pdb -o %s_BUF.pdb -ci buffer.pdb -nmol %s" %
         (pdbName, pdbName, countACID))

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

    file.write("BUF\t\t\t\t\t  %s\n" % (countACID))
topList.clear()

# GMX SOLVATE (ADD WATER) ######################################################

os.system("gmx solvate -cp %s_BUF.pdb -o %s_SOL.pdb -p topol.top" % (pdbName, pdbName))

# GMX GENION (ADD IONS) ########################################################

os.system("gmx grompp -f ions.mdp -c %s_SOL.pdb -p topol.top -o ions.tpr" % (pdbName))
os.system("gmx genion -s ions.tpr -o %s_ION.pdb -p topol.top -pname NA -nname CL -neutral << EOF\nSOL\nEOF" % pdbName)

# CREATE INDEX FILE ############################################################

group_AA  = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PRO', 'SER', 'THR', 'VAL', 'PHE', 'TYR', 'TRP']
group_BUF = ['BUF']
group_SOL = ['SOL']
group_ION = [' NA', ' CL']
group_SYS = group_AA  + group_BUF + group_SOL + group_ION
group_ISO = group_SOL + group_ION
# group_BASE   = ['ARG', 'LYS']
# group_ACID   = ['ASP', 'GLU']

protein2 = PDB()
protein2.loadpdb("%s_ION.pdb" % (pdbName))

protein2.writendx("index.ndx", "SYSTEM"    , group_SYS )
protein2.writendx("index.ndx", "PROTEIN"   , group_AA  )
protein2.writendx("index.ndx", "BUFFER"    , group_BUF )
protein2.writendx("index.ndx", "WATER"     , group_SOL )
protein2.writendx("index.ndx", "IONS"      , group_ION )
protein2.writendx("index.ndx", "WATER_IONS", group_ISO )
protein2.writendx("index.ndx", "ZZZ", ['ZZZ'] )

# ENERGY MINIMIZATION ##########################################################

# gmx grompp -f minim.mdp -c cardio_ions.pdb -p topol.top -o em.tpr
# gmx mdrun -s em.tpr -c em.pdb
