#!/usr/bin/python3

import os
from lib import *

# PARAMS #######################################################################

gromPath    = "/home/anton/GIT/phbuilder/grom"
modelFF     = "charmm36-mar2019"
modelWater  = "tip3p"

protein     = PDB("1cvo.pdb")

# CONSTANT-PH OPTIONS

pH          = 4.5           

# PATH AND FILE STUFF ##########################################################

# Copy files from grom to our working dir.
os.system("cp -r %s/* ." % gromPath)

# Backup previous builder.log file and index.ndx (if exists)
backupFile("builder.log")

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

# Create EOF string required for pdb2gmx to set the protonation state of 
# ASP and GLU to true (specify 1 for user input option.
xstr = "<< EOF"
for idx in range(0, countACID):
    xstr += "\n1"
xstr += "\nEOF"

print("           : Running gmx pdb2gmx to create %s_PR2.pdb..." % (pdbName))

# Generate topology and protonate (make neutral) all GLU and ASP:
os.system("gmx pdb2gmx -f %s_PR1.pdb -o %s_PR2.pdb -asp -glu -ignh -ff %s \
           -water %s >> builder.log 2>&1 %s" % (pdbName, pdbName, modelFF, 
           modelWater, xstr)) # Supress both stdout and stderr.

# GMX EDITCONF #################################################################

sim.protein.add.box(pdbName)        # Editconf

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

sim.protein.add.water(pdbName)      # Solvate

sim.protein.add.ions(pdbName)       # Genion

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

sim.generate.mdp("EM.mdp",  Type='EM',  dt=0.01,  nsteps=10000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate.mdp("NVT.mdp", Type='NVT', dt=0.002, nsteps=25000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate.mdp("NPT.mdp", Type='NPT', dt=0.002, nsteps=25000,  output=0,    tgroups=[['SYSTEM', 0.5, 300]])
sim.generate.mdp("MD.mdp",  Type='MD',  dt=0.002, nsteps=500000, output=1000, tgroups=[['SYSTEM', 0.5, 300]])

# CREATE CONSTANT-PH .DAT FILE #################################################

lambdaGen('%s_ION.pdb' % (pdbName), pH)

sim.write.run(pdbName, "/usr/local/gromacs", "/usr/local/gromacs_dev")
sim.write.reset(pdbName)
sim.write.jobscript(pdbName, "test", 36, 1)

# ENERGY MINIMIZATION ##########################################################

sim.energy.minimize(pdbName)
sim.energy.tcouple(pdbName)
sim.energy.pcouple(pdbName)
    