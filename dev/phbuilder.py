#!/usr/bin/python3

import os
from lib import *

# PARAMS #######################################################################

protein    = PDB("2khm.pdb", CHAIN='A')

gromPath   = "/home/anton/GIT/phbuilder/grom"   # relative path to grom dir
modelFF    = "charmm36-mar2019"
modelWater = "tip3p"

startProto = True                               # ASP and GLU will be neutral

# PATH AND FILE STUFF ##########################################################

# Copy files from grom to our working dir.
os.system("cp -r %s/* ." % gromPath)

# Backup previous builder.log file and index.ndx (if exists)
backupFile("builder.log")
backupFile("index.ndx")

# GMX PDB2GMX ##################################################################

pdbName = protein.fname(ext = 0)
# protein.resetResId()
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
# This piece of code is kind of a hoax but it works.
topList = []
with open("topol.top", "r") as file:
    for line in file.readlines():
        topList.append(line)

with open("topol.top", "w+") as file:
    try:
        for idx in range(0, len(topList)):
            file.write(topList[idx])

            # If we see that the next line is this:
            if topList[idx + 1] == "; Include water topology\n":
                # Then insert the buffer topology before that line:
                file.write("; Include buffer topology\n")
                file.write("#include \"buffer.itp\"\n\n")            

    except IndexError:
        pass

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

group_PROTEIN     = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'SER', 'THR',
                     'VAL', 'PHE', 'TYR', 'TRP']
group_BUFFER      = ['BUF']
group_WATER       = ['SOL']
group_IONS        = [' NA', ' CL']
group_NON_PROTEIN = group_BUFFER + group_WATER + group_IONS
group_SYSTEM      = group_PROTEIN + group_NON_PROTEIN

protein2 = PDB("%s_ION.pdb" % (pdbName))
protein2.writendx("index.ndx", "SYSTEM"     , group_SYSTEM      )
protein2.writendx("index.ndx", "PROTEIN"    , group_PROTEIN     )
protein2.writendx("index.ndx", "BUFFER"     , group_BUFFER      )
protein2.writendx("index.ndx", "WATER"      , group_WATER       )
protein2.writendx("index.ndx", "IONS"       , group_IONS        )
protein2.writendx("index.ndx", "NON_PROTEIN", group_NON_PROTEIN )

# CREATE .MDP FILES ############################################################

mdpGen("EM.mdp", Type='EM', dt=0.01, nsteps=10000, output=0,
       tgroups=[['SYSTEM', 0.5, 300]])

mdpGen("NVT.mdp", Type='NVT', dt=0.002, nsteps=25000, output=0,
       tgroups=[['SYSTEM', 0.5, 300]])

mdpGen("NPT.mdp", Type='NPT', dt=0.002, nsteps=25000, output=0,
       tgroups=[['SYSTEM', 0.5, 300]])

mdpGen("MD.mdp", Type='MD', dt=0.002, nsteps=500000, output=1000,
       tgroups=[['SYSTEM', 0.5, 300]])

lambdaGen('%s_ION.pdb' % (pdbName), 4.5)

# ENERGY MINIMIZATION ##########################################################

print("pHbuilder  : Running gmx grompp to create EM.tpr...")

os.system("gmx grompp -f EM.mdp -c %s_ION.pdb -p topol.top -n index.ndx \
           -o EM.tpr >> builder.log 2>&1" % (pdbName))

print("           : Running gmx mdrun (energy minimization) to create %s_EM.pdb..." % (pdbName))

os.system("gmx mdrun -s EM.tpr -o EM.trr -c %s_EM.pdb -g EM.log -e EM.edr \
           >> builder.log 2>&1" % (pdbName))

# TEMPERATURE COUPLING #########################################################

print("pHbuilder  : Running gmx grompp to create NVT.tpr...")

os.system("gmx grompp -f NVT.mdp -c %s_EM.pdb -p topol.top -n index.ndx \
           -o NVT.tpr -r %s_EM.pdb >> builder.log 2>&1" % (pdbName, pdbName))

print("           : Running gmx mdrun (temperature coupling) to create %s_NVT.pdb..." % (pdbName))

os.system("gmx mdrun -s NVT.tpr -o NVT.trr -c %s_NVT.pdb -g NVT.log -e NVT.edr \
           >> builder.log 2>&1" % (pdbName))

# PRESSURE COUPLING ############################################################

print("pHbuilder  : Running gmx grompp to create NPT.tpr...")

os.system("gmx grompp -f NPT.mdp -c %s_NVT.pdb -p topol.top -n index.ndx \
           -o NPT.tpr -r %s_NVT.pdb >> builder.log 2>&1" % (pdbName, pdbName))

print("           : Running gmx mdrun (pressure coupling) to create %s_NPT.pdb..." % (pdbName))

os.system("gmx mdrun -s NPT.tpr -o NPT.trr -c %s_NPT.pdb -g NPT.log -e NPT.edr \
           >> builder.log 2>&1" % (pdbName))

# PRODUCTION RUN ###############################################################

# print("pHbuilder  : Running gmx grompp to create MD.tpr...")

# os.system("gmx grompp -f MD.mdp -c %s_NPT.pdb -p topol.top -n index.ndx \
#            -o MD.tpr >> builder.log 2>&1" % (pdbName))

# os.system("gmx mdrun -v -s MD.tpr -o MD.trr -c %s_MD.pdb -g MD.log -e MD.edr" % (pdbName))

# OUR bug with 2khm.pdb is because of calibration, we use wrong ddlv values