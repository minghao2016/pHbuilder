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

    # Update index.ndx
    generate_index()

# MDP
################################################################################

def generate_mdp(Type, nsteps=25000, nstxout=0, posres=False):
    # HEAD
    if (Type not in ['EM', 'NVT', 'NPT', 'MD']):
        raise Exception("Unknown .mdp Type specified. Types are: IONS, EM, NVT, NPT, MD.")

    __update("generate_mdp", "Type={0}, nsteps={1}, nstxout={2}, posres={3}".format(Type, nsteps, nstxout, posres))

    file = open("{0}.mdp".format(Type), 'w')

    def addTitle(title):
        file.write("\n; {0}\n".format(title.upper()))

    def addParam(name, value, comment=''):
        if (comment == ''):
            file.write("{:20s} = {:13s}\n".format(name, str(value)))
        else:
            file.write("{:20s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    # POSITION RESTRAIN SECTION
    if Type in ['EM', 'MD']:
        if GET('d_constantpH') and GET('d_restrainpH') and posres:
            addTitle('Position restrain')
            addParam('define', '-DPOSRES -DPOSRES_BUF', 'Position restraints.')
        elif GET('d_constantpH') and GET('d_restrainpH'):
            addTitle('Position restrain')
            addParam('define', '-DPOSRES_BUF', 'Position restraints.')

    if (Type in ['NVT', 'NPT']): # position restrain temp and press coupling.
        addTitle('Position restrain')
        if GET('d_constantpH') and GET('d_restrainpH'):
            addParam('define', '-DPOSRES -DPOSRES_BUF', 'Position restraints.')
        else:
            addParam('define', '-DPOSRES', 'Position restraints.')

    # RUN CONTROL
    addTitle("Run control")

    if (Type == 'EM'): # emtol hardcored, pretty typical for normal MD.
        dt = 0.01
        addParam('integrator', 'steep', 'Use steep for EM.')
        addParam('emtol', 1000, 'Stop when max force < 1000 kJ/mol/nm.')
        addParam('emstep', dt, 'Time step (ps).')

    if (Type in ['NVT', 'NPT', 'MD']):
        dt = 0.002
        addParam('integrator', 'md')
        addParam('dt', dt, 'Time step (ps).')

    addParam('nsteps', nsteps, '%.1f ns.' % ((dt * nsteps)/1000.0))

    # We restrain the COM to prevent protein from coming too close to the BUFs.
    if Type == 'MD' and GET('d_restrainpH'):
        addParam('comm-mode', 'Linear', 'Remove center of mass translation.')
        addParam('comm-grps', 'Protein Non-Protein')

    # OUTPUT CONTROL
    addTitle("Output control")
    addParam('nstxout-compressed', nstxout, 'Write frame every %.3f ps.' % (dt * nstxout))

    # NEIGHBOUR SEARCHING PARAMETERS
    addTitle("Neighbour searching")
    addParam('cutoff-scheme', 'Verlet', 'Related params are inferred by Gromacs.')

    # BONDED
    if (Type in ['NVT', 'NPT', 'MD']):
        addTitle("Bond parameters")
        addParam('constraints', 'h-bonds', 'Constrain H-bond vibrations.')
        addParam('constraint_algorithm', 'lincs', 'Holonomic constraints.')
        addParam('lincs_iter', 1, 'Related to accuracy of LINCS.')
        addParam('lincs_order', 4, 'Related to accuracy of LINCS.')

    # ELECTROSTATICS
    addTitle("Electrostatics")
    addParam('coulombtype', 'PME', 'Use Particle Mesh Ewald.')
    
    if (GET('d_modelFF')[0:5].lower() == "charm"): # if we use a CHARMM force field...
        addParam('rcoulomb', 1.2, 'Berk: CHARMM was calibrated for 1.2 nm.')
        addParam('fourierspacing', 0.14, 'Berk: set this to 0.14 for CHARMM.')
    else: # Default for force fields:
        addParam('rcoulomb', 1.0, 'Coulomb cut-off (nm).')

    # VAN DER WAALS
    addTitle("Van der Waals")
    addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')

    if (GET('d_modelFF')[0:5].lower() == "charm"): # if we use a CHARMM force field...
        addParam('rvdw', 1.2, 'Berk: CHARMM was calibrated for 1.2 nm.')
        addParam('vdw-modifier', 'force-switch', 'Berk: specific for CHARMM.')
        addParam('rvdw-switch', 1.0, 'Berk: specific for CHARMM.')
    else: # Default for force fields:
        addParam('rvdw', 1.0, 'Van der Waals cut-off (nm).')

    # TEMPERATURE COUPLING
    tgroups = [['SYSTEM', 0.5, 300]]

    if (Type in ['NVT', 'NPT', 'MD']):
        addTitle("Temperature coupling")
        addParam('tcoupl', 'v-rescale')

        string1 = ""; string2 = ""; string3 = ""
        for group in tgroups:
            string1 += group[0]      + ' '
            string2 += str(group[1]) + ' '
            string3 += str(group[2]) + ' '

        addParam('tc-grps', string1)
        addParam('tau-t', string2, 'Berk: change from 0.1 to 0.5.')
        addParam('ref-t', string3, 'Reference temp. (K) (for each group).')

    # PRESSURE COUPLING
    if (Type in ['NPT', 'MD']):
        addTitle('Pressure coupling')
        
        if (Type == 'NPT'):
            addParam('pcoupl', 'Berendsen', 'Use Berendsen for NPT.')
        else:
            addParam('pcoupl', 'Parrinello-Rahman')
        
        addParam('pcoupltype', 'isotropic', 'Uniform scaling of box.')
        addParam('tau_p', 5.0, 'Berk: better to change from 2.0 to 5.0.')
        addParam('ref_p', 1.0, 'Reference pressure (bar).')
        addParam('compressibility', 4.5e-5, 'Isothermal compressbility of water.')
        addParam('refcoord_scaling', 'all', 'Required with position restraints.')

    # PERIODIC BOUNDARY CONDITIONS
    addTitle("Periodic boundary condition")
    addParam('pbc', 'xyz', 'To keep molecule(s) in box.')

    # GENERATE VELOCITIES FOR STARTUP
    if (Type == 'NVT'):
        addTitle('Generate velocities for startup')
        addParam('gen_vel', 'yes')

    file.close()

    # PUT RELEVANT PARAMETERS IN UNIVERSE
    if (Type == 'MD'):
        UPD('d_dt', dt)
        UPD('d_nsteps', nsteps)
        UPD('d_nstxout', nstxout)

# TOPOLOGY/SIMULATION BUILDER STUFF
################################################################################

def add_mol_to_topol(itpfname, comment, molname=None, molcount=None):
    # Get the contents of current topol.top.
    topList = []
    with open("topol.top") as file:
        for line in file.readlines():
            topList.append(line)
    
    # Add the .itp file (line saying: #include "blabla.itp")
    with open("topol.top", 'w') as file:
        try:
            for line in range(0, len(topList)):
                file.write(topList[line])
                
                if "[ system ]\n" == topList[line + 1]:
                    file.write("; {0}\n".format(comment))
                    file.write("#include \"{0}\"\n\n".format(itpfname))

        except IndexError:
            pass

    # if molcount not present, add it, otherwise do nothing.
        if molname != None and molcount != None and molname not in topList[-1]:
            file.write("{0}\t\t\t{1}\n".format(molname, molcount))

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
  
    with open('topol.top', 'w') as file:
        file.write("; Include forcefield parameters\n")
        file.write("#include \"{0}.ff/forcefield.itp\"\n\n".format(GET('d_modelFF')))

        file.write("; Include protein topology\n")
        for letter in GET('d_chain'):
            file.write("#include \"topol_Protein_chain_{0}.itp\"\n".format(letter))
        file.write('\n')

        file.write('[ system ]\n')
        # file.write('; name\n')
        file.write('{0}\n\n'.format(GET('d_pdbName')))

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

    # Rebuild topology
    __rebuild_topol()

    # Update d_residues.
    loadpdb("{0}_PR2.pdb".format(GET('d_pdbName')))

    # Update index.ndx
    generate_index()

def protein_add_box(d_boxMargin, d_boxType='cubic'):
    __update("protein_add_box", "running gmx editconf...")

    os.system("gmx editconf -f {0}_PR2.pdb -o {0}_BOX.pdb -d {1} -bt {2} >> builder.log 2>&1".format(GET('d_pdbName'), d_boxMargin, d_boxType))

    SET(d_boxMargin)
    SET(d_boxType)

def protein_add_buffer(bufpdbName, bufitpname, minSep=1.5):
    if not (GET('d_constantpH') and GET('d_restrainpH')):
        __update("protein_add_buffer", "skipping this step...")
        return
    
    __update("protein_add_buffer", "adding buffer molecules...")

    def extractMinimum():                  # Extract the required value.
        def Float(fileName, line, col):
            for x, y in enumerate(open(fileName)):
                if (x == line - 1):
                    return float(y.split()[col-1])
                                                
        return Float("mindist.xvg", 25, 2) # Position of mindist in .xvg file.

    def clean():
        os.system("rm -f \\#mindist.xvg.*\\# \\#{0}_BUF.pdb.*\\# mindist.xvg".format(GET('d_pdbName')))

    countACID = countRes('ASP') + countRes('GLU')

    attempts = 0
    while (True):           # Randomly add the buffer molecules
        os.system("gmx insert-molecules -f {0}_BOX.pdb -o {0}_BUF.pdb -ci {1} -nmol {2} >> builder.log 2>&1".format(GET('d_pdbName'), bufpdbName, countACID))
                            # Determine the minimum distance between the protein and the buffers
        os.system("gmx mindist -f {0}_BUF.pdb -s {0}_BUF.pdb >> builder.log 2>&1 << EOF\n1\n13\nEOF".format(GET('d_pdbName')))

        attempts += 1
        if (extractMinimum() >= minSep):
            break

        if (attempts > 100):
            clean()
            raise Exception("Maximum number of buffer insertion attempts exceeded (100). Try decreasing minSep or increasing boxSizeMargin.")

    __update("protein_add_buffer", "placing buffers took %s attempt(s) (mindist = %s)..." % (attempts, extractMinimum()))
    
    clean()

    __update("protein_add_buffer", "updating topology...")
    
    # Add buffer topology to topol.top.
    add_mol_to_topol(bufitpname, "Include buffer topology", 'BUF', countACID)

    # Update d_residues.
    loadpdb("{0}_BUF.pdb".format(GET('d_pdbName'))) 

    # Update index.ndx
    generate_index()

def protein_add_water():
    __update("protein_add_water", "running gmx solvate...")
    
    # Add water molecules to .pdb and update [ molecules ] in .top:
    if (GET('d_constantpH') and GET('d_restrainpH')):
        os.system("gmx solvate -cp {0}_BUF.pdb -o {0}_SOL.pdb -p topol.top >> builder.log 2>&1".format(GET('d_pdbName')))
    else:
        os.system("gmx solvate -cp {0}_BOX.pdb -o {0}_SOL.pdb -p topol.top >> builder.log 2>&1".format(GET('d_pdbName')))

    # Update topol.top
    add_mol_to_topol("{0}.ff/{1}.itp".format(GET('d_modelFF'), GET('d_modelWater')), "Include water topology", "SOL", countRes('SOL'))

    loadpdb("{0}_SOL.pdb".format(GET('d_pdbName'))) # Update d_residues.

    # Update index.ndx
    generate_index()

def protein_add_ions():
    __update("protein_add_ions", "running gmx grompp and genion to add ions...")

    # Generate IONS.mdp (just a dummy required).
    os.system('touch IONS.mdp')

    # Add ion topology to topol.top.
    add_mol_to_topol("{0}.ff/ions.itp".format(GET('d_modelFF')), "Include ion topology")

    # Run gmx genion.
    os.system("gmx grompp -f IONS.mdp -c {0}_SOL.pdb -p topol.top -o IONS.tpr >> builder.log 2>&1".format(GET('d_pdbName')))
    os.system("gmx genion -s IONS.tpr -o {0}_ION.pdb -p topol.top -pname NA -nname CL -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF".format(GET('d_pdbName')))

    # Update d_residues.
    loadpdb("{0}_ION.pdb".format(GET('d_pdbName')))

    # Update index.ndx
    generate_index()

# INDEX, JOBSCRIPT, RUN, RESET
################################################################################
def generate_index():
    #__update("generate_index", "running gmx make_ndx to create index.ndx...")
    os.system("gmx make_ndx -f {0}_ION.pdb -o index.ndx >> builder.log 2>&1 << EOF\nq\nEOF".format(GET('d_pdbName')))

def write_run(gmxPath="/usr/local/gromacs", options=""):
    __update("write_run", "gmxPath={0}, additional options= {1}".format(gmxPath, options))
    
    with open("run.sh", 'w') as file:
        file.write("#!/bin/bash\n\n")

        file.write("# Gromacs version to use:\n")
        file.write("source {0}/bin/GMXRC\n\n".format(gmxPath))

        file.write("gmx grompp -f MD.mdp -c {0}_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r {0}_NPT.pdb\n".format(GET('d_pdbName')))
        file.write("gmx mdrun -v -s MD.tpr -x MD.xtc -c {0}_MD.pdb -g MD.log -e MD.edr {1}\n".format(GET('d_pdbName'), options))

    os.system("chmod +x run.sh")

def write_jobscript(jobName, jobTime, nodes, ntasks, queue):
    __update("write_jobscript", "jobName={0}, jobTime={1}(hrs), nodes={2}, ntasks={3}, queue={4}...".format(jobName, jobTime, nodes, ntasks, queue))

    file = open("jobscript.sh", 'w')
    
    def writeHead(param, value):
        file.write("#SBATCH --%s=%s\n" % (param, value))

    def moduleLoad(value):
        file.write("module load {0}\n".format(value))

    file.write("#!/bin/bash\n")

    writeHead("time", "%d-%.2d:00:00" % (int(jobTime / 24), jobTime % 24))
    writeHead("nodes", nodes)
    writeHead("ntasks", ntasks)
    writeHead("partition", queue)
    writeHead("job-name", jobName)
    writeHead("mail-user", "anton.jansen@scilifelab.se")
    writeHead("mail-type", "ALL")
    writeHead("C gpu --gres", "gpu:1")
    
    file.write('\n')

    moduleLoad("cmake/latest")
    moduleLoad("cuda/10.2")

    file.write('\n')

    if GET('d_constantpH'):
        file.write("\n# compile our custom Gromacs version on cluster backend node\n")
        file.write("mkdir build\n")
        file.write("cd build\n")
        file.write("cmake ~/gromacs-constantph -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=${PWD}/.. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA\n")
        file.write("make -j 12\n")
        file.write("make install -j 12\n")
        file.write("cd ..\n")
        file.write("rm -r build\n")
        file.write("source ${PWD}/bin/GMXRC\n\n")
    else:
        file.write("module load gromacs/2021.1\n\n")

    file.write("gmx grompp -f MD.mdp -c {0}_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r {0}_NPT.pdb\n".format(GET('d_pdbName')))
    
    if GET('d_constantpH'):
        file.write("gmx mdrun -s MD.tpr -x MD.xtc -c {0}_MD.pdb -g MD.log -e MD.edr -pme cpu\n".format(GET('d_pdbName')))
    else:
        file.write("gmx mdrun -s MD.tpr -x MD.xtc -c {0}_MD.pdb -g MD.log -e MD.edr\n".format(GET('d_pdbName')))

def write_reset():
    __update("write_reset", "writing reset.sh...")

    with open("reset.sh", "w+") as file:
        file.write("#!/bin/bash\n\n")
        
        file.write("if [ -f \"%s_MD.pdb\" ]\nthen\n" % GET('d_pdbName'))
        file.write("\tread -p \"Warning: simulation has finished. Proceed? (y)\" var\n")
        file.write("else\n")
        file.write("\trm -rf \\_\\_py* charmm*\n")
        file.write("\trm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg\n")
        file.write("\trm -f \\#*\\#\n")
        file.write("\trm -f buffer.pdb %s_*.pdb\n" % GET('d_pdbName'))
        file.write("\trm -f run.sh reset.sh jobscript.sh universe\n")                   
        file.write("fi\n\n")

        file.write("if [ \"${var}\" = \"y\" ]\nthen\n")
        file.write("\trm -rf \\_\\_py* charmm*\n")
        file.write("\trm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg\n")
        file.write("\trm -f \\#*\\#\n")
        file.write("\trm -f buffer.pdb %s_*.pdb\n" % GET('d_pdbName'))
        file.write("\trm -f run.sh reset.sh jobscript.sh universe\n")            
        file.write("fi\n\n")

    os.system("chmod +x reset.sh")    

# ENERGY MINIMIZATION/COUPLING
################################################################################

def energy_minimize():
    generate_mdp('EM') # Create .mdp file for EM run.

    __update("energy_minimize", "running gmx grompp and mdrun for energy minimization...")

    os.system("gmx grompp -f EM.mdp -c {0}_ION.pdb -p topol.top -n index.ndx -o EM.tpr -r {0}_ION.pdb >> builder.log 2>&1".format(GET('d_pdbName')))
    os.system("gmx mdrun -s EM.tpr -o EM.trr -c {0}_EM.pdb -g EM.log -e EM.edr >> builder.log 2>&1".format(GET('d_pdbName')))

def energy_tcouple():
    generate_mdp('NVT') # Create .mdp file for temperature coupling run.
    
    __update("energy_tcouple", "running gmx grompp and mdrun for temperature coupling...")

    os.system("gmx grompp -f NVT.mdp -c {0}_EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r {0}_EM.pdb >> builder.log 2>&1".format(GET('d_pdbName')))
    os.system("gmx mdrun -s NVT.tpr -o NVT.trr -c {0}_NVT.pdb -g NVT.log -e NVT.edr >> builder.log 2>&1".format(GET('d_pdbName')))

def energy_pcouple(skip=False):
    generate_mdp('NPT') # Create .mdp file for pressure coupling run.
    __update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")

    if (skip): # This can save some time in some cases.
        __update("energy_pcouple", "skipping pressure coupling (only copying .pdb)...")
    
        os.system("cp {0}_NVT.pdb {0}_NPT.pdb".format(GET('d_pdbName')))
    else:
        __update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")

        os.system("gmx grompp -f NPT.mdp -c {0}_NVT.pdb -p topol.top -n index.ndx -o NPT.tpr -r {0}_NVT.pdb >> builder.log 2>&1".format(GET('d_pdbName')))
        os.system("gmx mdrun -s NPT.tpr -o NPT.trr -c {0}_NPT.pdb -g NPT.log -e NPT.edr >> builder.log 2>&1".format(GET('d_pdbName')))    

# GENERATE CONSTANT-PH DATA
################################################################################

def generate_phdata(pH, lambdaM, nstOut, barrierE):
    pass

# MAIN #########################################################################

UPD('d_constantpH', True)
UPD('d_restrainpH', True)

processpdb('3lr2.pdb')
write_reset()

generate_topology("charmm36-mar2019", "tip3p")
protein_add_box(d_boxMargin=2.0)
protein_add_buffer("/home/anton/GIT/phbuilder/grom/buffer.pdb", "/home/anton/GIT/phbuilder/grom/buffer.itp", 1.5)
protein_add_water()
protein_add_ions()

generate_index()

write_run()
write_jobscript('test', 48, 1, 32, 'lindahl')

generate_mdp('MD', nsteps=500000, nstxout=10000)

energy_minimize()
energy_tcouple()
energy_pcouple(skip=True)

# PRINT()
