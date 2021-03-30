#!/bin/python3

import shelve       # This module allows retrieving of variable names.
import varname      # https://github.com/pwwang/python-varname
import os

def PRINT(): # Displays all variables (name, data, type) stored in the universe.
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in shelf:
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}]  ({3})".format(item, shelf[item][0], shelf[item][-1], type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1}  ({2})".format(item, shelf[item], type(shelf[item]), arg=longest))

def SET(var): # Add/update variable to universe.
    with shelve.open('universe') as shelf:
        shelf[varname.nameof(var, frame=2)] = var
        # Not sure how stable this is, thread carefully!

def UPD(varName, value): # Update variable to universe.
    with shelve.open('universe') as shelf:
        shelf[varName] = value

def GET(var): # Retrieve variable from universe.
    with shelve.open('universe') as shelf:
        try:
            return shelf[var]
        except KeyError:
            data = eval(input("GET couldn't retrieve var \"{0}\" from {1}. Enter manually: ".format(var, 'universe')))
            shelf[var] = data
            print("SET {0} = {1}  ({2})".format(var, shelf[var], type(shelf[var])))
            return shelf[var]

def __update(tool, message): # For formating user messages.
    print("{:22s} : {:s}".format(tool, message))

# PDB STUFF
################################################################################

class Residue: # Stores a single residue's data.
    def __init__(self, atoms, ali, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms      # list      holds atom names
        self.d_ali     = ali        # list      holds alternative loc. ind.
        self.d_resname = resname    # string    holds residue name
        self.d_chain   = chain      # string    holds chain name (A, B)
        self.d_resid   = resid      # int       holds residue number
        self.d_x       = x          # list      holds x-coordiantes
        self.d_y       = y          # list      holds y-coordinates
        self.d_z       = z          # list      holds z-coordinatees

# Load a .pdb file into the universe.
def loadpdb(d_fname, d_model=1, d_ALI='A', d_chain=["all"]):
    d_residues = []
    atomLines  = []
    with open(d_fname) as file:     # Read .pdb line-by-line.
        read = True                 # True if no specific MODEL specified.

        for line in file.readlines():
            # This is to make sure we only import the specified MODEL.
            if ((line[0:6]) == "MODEL "):
                if ("MODEL {:8d}".format(d_model) in line):
                    read = True
                else:
                    read = False
            
            # Get title.
            if ((line[0:6]) == "TITLE "):
                d_title = line[7:80].rstrip(); SET(d_title)

            # Get periodic box information (if any).
            if ((line[0:6]) == "CRYST1"):
                d_box   = line[7:80].rstrip(); SET(d_box)
                
            # if our line specifies an ATOM,
            if (line[0:6] == "ATOM  "):                 
                # and we are currently reading the correct MODEL,
                if (read == True):
                    # and our line contains the correct specified alternate 
                    # location specifier (if any at all),
                    if (line[16:17] in [d_ALI, " "]):
                        # and we want all chains,
                        if (d_chain == ["all"]):
                            # then load the atom.
                            atomLines.append(line)
                        # or if we want a selection of chains,
                        elif (line[21:22] in d_chain):
                            # load that selection.
                            atomLines.append(line)

    # Add one line of padding to prevent IndexError
    atomLines.append("0000000000000000000000000000000000000000000000000000")

    # Loop through lines and create list of Residue objects.
    atoms = []; ali = []; xCoord = []; yCoord = []; zCoord = []
    
    for idx in range(0, len(atomLines) - 1):
        atoms.append(atomLines[idx][12:16])
        ali.append(atomLines[idx][16:17])
        xCoord.append(float(atomLines[idx][30:38]))
        yCoord.append(float(atomLines[idx][38:46]))
        zCoord.append(float(atomLines[idx][46:54]))

        # If the resid of the next line is different, we are at end
        if (atomLines[idx][22:26] != atomLines[idx + 1][22:26]):
            # put the data in a Residue object and append to d_residues:
            d_residues.append(
                Residue
                (
                    atoms, 
                    ali,
                    atomLines[idx][17:20], 
                    atomLines[idx][21:22], 
                    int(atomLines[idx][22:26]), 
                    xCoord, 
                    yCoord,
                    zCoord
                ))
            # Reset.
            atoms = []; ali = []; xCoord = []; yCoord = []; zCoord = []

    SET(d_residues)

# Write (modified) .pdb file.
def __writepdb(name):            
    with open(name, 'w') as file:
        file.write("TITLE {0}\n".format(GET('d_title')))
        file.write("CRYST1{0}\n".format(GET('d_box')))
        file.write("MODEL {:8d}\n".format(GET('d_model')))

        num = 1                         # Write actual residues.
        for residue in GET('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', num, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                num += 1                # .pdb format string.

        file.write("TER\nENDMDL\n")     # Write EOF information.

# Inspect a specific residue.
def inspectRes(resid, chain='A'):
    for residue in GET('d_residues'):
        if (residue.d_resid == resid and residue.d_chain == chain):
            print("resname = {0}, resid = {1}, chain = {2}".format(residue.d_resname, residue.d_resid, residue.d_chain))

            for idx in range(0, len(residue.d_atoms)):
                print("{:^4s}{:1s}{:8.3f} {:8.3f} {:8.3f}".format(residue.d_atoms[idx], residue.d_ali[idx], residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))

            return

    raise Exception("Cannot find residue with resid={0} and/or chain={1}".format(resid, chain))

# Count number of residues with a specific residue name.
def countRes(resname):
    count = 0
    for residue in GET('d_residues'):
        if (residue.d_resname == resname):
            count += 1
    
    return count

# Reset the original residue numbering.
def __resetResId():
    num = 1
    for residue in GET('d_residues'):
        residue.d_resid = num
        num += 1

def __genMissingChainIdentifiers():
    import string
    # This function is based on the fact that the resid resets for every
    # separate chain. If that does not happen, this won't work.
    capIdx = 0
    caps   = string.ascii_uppercase

    for idx in range(0, len(GET('d_residues'))):
        try:
            current = GET('d_residues')[idx]
            Next    = GET('d_residues')[idx + 1]

            GET('d_residues')[idx].d_chain = caps[capIdx]

            # If the resid resets, we know we're at the end of a chain.
            if (current.d_resid + 1 != Next.d_resid):
                capIdx += 1
        
        except IndexError:
            GET('d_residues')[idx].d_chain = caps[capIdx]    

def processpdb(d_fname, d_model=1, d_ALI='A', d_chain=["all"], resetResId=False, genMissingChainsIDs=False):
    SET(d_fname)
    SET(d_model)
    SET(d_ALI)
    
    loadpdb(d_fname, d_model, d_ALI, d_chain)       # Load.

    d_pdbName = d_fname[0:len(d_fname)-4]           # Set internal name.
    SET(d_pdbName)

    d_chain = []                                    # Update d_chain from ["all"] 
    for residue in GET('d_residues'):               # to a list of actual chains used.
        d_chain.append(residue.d_chain)
    d_chain = list(set(d_chain))
    d_chain.sort()
    SET(d_chain)

    if (genMissingChainsIDs):                       # Optional stuff.
        __genMissingChainIdentifiers()

    if (resetResId):
        __resetResId()
                                                    # Write processed file.
    __writepdb("{0}_PR1.pdb".format(GET('d_pdbName')))

# TOPOLOGY/SIMULATION BUILDER STUFF
################################################################################

def __rebuild_topol():
    # If we have only one chain, gromacs will put everything in topol.top.
    # If we have more than one chain, gromacs will do it for us.
    if (len(GET('d_chain')) <= 1):
        readingProtein = False
        
        file = open("topol_Protein_chain_A.itp", 'w')

        for line in open("topol.top").readlines():
            if (not readingProtein and line == "[ moleculetype ]\n"):
                readingProtein = True
        
            if (readingProtein and line == "; Include water topology\n"):
                readingProtein = False

            if (readingProtein):
                file.write(line)
   
        file.close()
        # file2.close()
  
    with open('topol.top', 'w') as file:
        file.write("; Include forcefield parameters\n")
        file.write("#include \"{0}.ff/forcefield.itp\"\n\n".format(GET('d_modelFF')))

        file.write("; Include protein topology\n")
        for letter in GET('d_chain'):
            file.write("#include topol_Protein_chain_{0}.itp\n".format(letter))
        file.write('\n')

        file.write('[ system ]\n')
        file.write('; name\n')
        file.write('TITLE\n\n')

        file.write('[ molecules ]\n')
        file.write('; Compounts \t\t #mols\n')
        for letter in GET('d_chain'):
            file.write("Protein_chain_{0}\t\t1\n".format(letter))

def generate_topology(d_modelFF, d_modelWater, d_neutralTermini=False):
    SET(d_modelFF)
    SET(d_modelWater)
    SET(d_neutralTermini)
    
    countACID = countRes("ASP") + countRes("GLU")

    # If constant-pH is on,
    if GET('d_constantpH'):
        __update("protein_add_forcefield", 'constant-pH is turned on...')
        
        # and we have at least one protonatable reside,
        if (countACID > 0):
            __update("protein_add_forcefield", "detected {0} acidic residue(s):".format(countACID))

            count = 1
            for residue in GET('d_residues'):
                if (residue.d_resname in ['ASP', 'GLU']):
                    __update("protein_add_forcefield", "{:3s} {:<4d}".format(residue.d_resname, count))
                count += 1
            __update("protein_add_forcefield", "(setting protonation state to true for all of these)")

        else:
            __update("protein_add_forcefield", "no acidic residues detected, turning off constant-pH...")
            UPD('d_constantpH', False)
            UPD('d_restrainpH', False)

    else:
        __update("protein_add_forcefield", 'constant-pH is turned off...')

    __update("protein_add_forcefield", "using the %s force field with the %s water model..." % (d_modelFF, d_modelWater))

    # Count how many different chains we have
    listChains = []
    # atomCount  = 1
    for residue in GET('d_residues'):
        # throw if a residue/atom does not have any chain identifier
        # if residue.d_chain == ' ':
            # string = "residue {} (atom {}) does not belong to any chain (does not have a chain-identifier)!".format(residue.d_resid, atomCount)
            # raise Exception(string)
        # for _ in residue.d_atoms:
            # atomCount += 1

        if residue.d_chain not in listChains:
            listChains.append(residue.d_chain)

    # Create EOF string to circumvent pd2gmx prompting for user input.
    # The for-loop sets the protonation state of all GLUs and ASPs to 1 
    # (True) if d_constantpH is True. The line below sets the termini
    # to 0: NH3+, 0: COO- (charged, gmx default) if d_neutralTermini is False,
    # and to 1: NH2, 1: COOH (neutral) if d_neutralTermini is True.
    xstr = "<< EOF"
    for _ in range(0, countACID):
        xstr += "\n%d" % GET('d_constantpH')
    for _ in range(0, len(listChains)):
        xstr += "\n%d\n%d" % (d_neutralTermini, d_neutralTermini)
    xstr += "\nEOF"

    if (d_neutralTermini):
        __update("protein_add_forcefield", "setting termini to neutral (NH2 and COOH)...")
    else:
        __update("protein_add_forcefield", "setting termini to charged (NH3+ and COO-) (gmx default)...")
    
    __update("protein_add_forcefield", "running pdb2gmx to create topol.top...")

    # Generate topology and protonate (make neutral) all GLU and ASP:
    if (d_neutralTermini):
        os.system("gmx pdb2gmx -f {0}_PR1.pdb -o {0}_PR2.pdb -asp -glu -ignh -ff {1} -water {2} -ter >> builder.log 2>&1 {3}".format(GET('d_pdbName'), d_modelFF, d_modelWater, xstr))
    else: # Still a bug here if we have multiple chains, so have to specify the termini multiple times
        os.system("gmx pdb2gmx -f {0}_PR1.pdb -o {0}_PR2.pdb -asp -glu -ignh -ff {1} -water {2} >> builder.log 2>&1 {3}".format(GET('d_pdbName'), d_modelFF, d_modelWater, xstr))

    # Update internal d_residues.
    loadpdb("%s_PR2.pdb" % GET('d_pdbName'))

    # Cleanup
    os.system('rm -f posre*.itp')

    # Rebuild topology
    __rebuild_topol()

def protein_add_box(d_boxMargin, d_boxType='cubic'):
    __update("protein_add_box", "running gmx editconf...")

    os.system("gmx editconf -f {0}_PR2.pdb -o {0}_BOX.pdb -d {1} -bt {2} >> builder.log 2>&1".format(GET('d_pdbName'), d_boxMargin, d_boxType))

    SET(d_boxMargin)
    SET(d_boxType)

    loadpdb("{0}_BOX.pdb".format(GET('d_pdbName'))) # Update internal d_residues.

# MAIN #########################################################################

processpdb('GluCl.pdb')
UPD('d_constantpH', True)
generate_topology("charmm36-mar2019", "tip3p")
# protein_add_box(d_boxMargin=1.0)

PRINT()




























# def __load_topol():
#     d_topol = []

#     with open("topol.top") as file:
#         for line in file.readlines():
#             d_topol.append(line)

#     SET(d_topol)

# def add_mol_to_topol(add_before_line_containing, fname, comment, molname, molcount):
#     __load_topol()
    
#     with open("topol.top", 'w+') as file:
#         try:
#             for idx in range(0, len(GET('d_topol'))):
#                 file.write(GET('d_topol')[idx])

#                 if add_before_line_containing in GET('d_topol')[idx + 1]:
#                     file.write("; {0}\n".format(comment))
#                     file.write("#include \"{0}\"\n\n".format(fname))

#         except IndexError:
#             pass

#         file.write("{0}\t\t\t{1}\n".format(molname, molcount))

#     __load_topol()
