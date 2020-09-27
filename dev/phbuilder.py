#!/usr/bin/python3

import os
from lib import *

# PARAMS #######################################################################

protein    = PDB("1cvo.pdb")            # Additional options when importing

gromPath   = "/home/anton/GIT/phbuilder/grom"   # relative path to grom dir
modelFF    = "charmm36-mar2019"
modelWater = "tip3p"

startProto = True                               # ASP and GLU will be neutral



# PATH AND FILE STUFF ##########################################################

# Copy files from grom to our working dir.
os.system("cp -r %s/* ." % gromPath)

# Backup previous builder.log file (if it exists)
backupFile("builder.log")

# GMX PDB2GMX ##################################################################

pdbName = protein.fname(ext = 0)
protein.writepdb("%s_PR1.pdb" % (pdbName))

countASP  = protein.countRes("ASP")
countGLU  = protein.countRes("GLU")
countACID = countASP + countGLU

# Print how many acidic residues were found
print("pHbuilder  : Detected %s acidic residues (%s ASP and %s GLU)..." 
      % (countACID, countASP, countGLU))

print("           : Setting initial protonation state to %s" % (startProto))

# Create EOF string required for pdb2gmx
xstr = "<< EOF"
for idx in range(0, countACID):
    xstr += "\n%s" % (int(startProto))
xstr += "\nEOF"

print("           : Running gmx pdb2gmx to create %s_PR2.pdb..." % (pdbName))

# Generate topology and protonate (make neutral) all GLU and ASP:
os.system("gmx pdb2gmx -f %s_PR1.pdb -o %s_PR2.pdb -asp -glu -ignh -ff %s \
           -water %s >> builder.log 2>&1 %s" % (pdbName, pdbName, modelFF, 
           modelWater, xstr)) # Supress both stdout and stderr.

# GMX EDITCONF #################################################################

print("pHbuilder  : Running gmx editconf to create %s_BOX.pdb..." % (pdbName))

os.system("gmx editconf -quiet -f %s_PR2.pdb -o %s_BOX.pdb -c -d 1.0 -bt cubic \
           >> builder.log 2>&1" % (pdbName, pdbName)) #Supress only stdout.

# GMX INSERT-MOLECULES (ADD BUFFER) ############################################

print("pHbuilder  : Running gmx insert-molecules to create %s_BUF.pdb..." % (pdbName))

os.system("gmx insert-molecules -quiet -f %s_BOX.pdb -o %s_BUF.pdb -ci buffer.pdb \
           -nmol %s >> builder.log 2>&1" % (pdbName, pdbName, countACID))

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

print("pHbuilder  : Running gmx solvate to create %s_BUF.pdb..." % (pdbName))

os.system("gmx solvate -cp %s_BUF.pdb -o %s_SOL.pdb -p topol.top \
           >> builder.log 2>&1" % (pdbName, pdbName))

# GMX GENION (ADD IONS) ########################################################

print("pHbuilder  : Running gmx grompp to create IONS.tpr...")

os.system("gmx grompp -f IONS.mdp -c %s_SOL.pdb -p topol.top -o IONS.tpr \
           >> builder.log 2>&1" % (pdbName))

print("           : Running gmx genion to create %s_ION.pdb..." % (pdbName))

os.system("gmx genion -s IONS.tpr -o %s_ION.pdb -p topol.top -pname NA \
           -nname CL -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF" % pdbName)

# CREATE INDEX FILE ############################################################

print("pHbuilder  : Creating index file...")

group_AA  = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
    'LEU', 'LYS', 'MET', 'PRO', 'SER', 'THR', 'VAL', 'PHE', 'TYR', 'TRP']
group_BUF = ['BUF']
group_SOL = ['SOL']
group_ION = [' NA', ' CL']
group_SYS = group_AA  + group_BUF + group_SOL + group_ION
group_ISO = group_SOL + group_ION

protein2 = PDB("%s_ION.pdb" % (pdbName))
protein2.writendx("index.ndx", "SYSTEM"    , group_SYS )
protein2.writendx("index.ndx", "PROTEIN"   , group_AA  )
protein2.writendx("index.ndx", "BUFFER"    , group_BUF )
protein2.writendx("index.ndx", "WATER"     , group_SOL )
protein2.writendx("index.ndx", "IONS"      , group_ION )
protein2.writendx("index.ndx", "WATER_IONS", group_ISO )

# # ENERGY MINIMIZATION ##########################################################

mdpGen("EM.mdp")

print("pHbuilder  : Running gmx grompp to create EM.tpr...")

os.system("gmx grompp -f EM.mdp -c %s_ION.pdb -p topol.top -n index.ndx -o EM.tpr  >> builder.log 2>&1" % (pdbName))

print("           : Running gmx mdrun (energy minimization) to create %s_EM.pdb..." % (pdbName))

os.system("gmx mdrun -s EM.tpr -c %s_EM.pdb >> builder.log 2>&1" % (pdbName))
