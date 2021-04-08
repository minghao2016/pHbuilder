#!/bin/python3

import shelve
import os

# Set/update variable to universe.
def SET(varName, value):
    with shelve.open('universe') as shelf:
        shelf[varName] = value

# Retrieve variable from universe.
def GET(varName):
    with shelve.open('universe') as shelf:
        try:
            return shelf[varName]
        except KeyError:
            data = eval(input("GET couldn't retrieve var \"{0}\" from {1}. Enter manually: ".format(varName, 'universe')))
            shelf[varName] = data
            print("SET {0} = {1}  {2}".format(varName, shelf[varName], type(shelf[varName])))
            return shelf[varName]

# Display all variables (name, data, type) stored in the universe.
def PRINT():
    with shelve.open('universe') as shelf:
        # Determine longest valueName for formatting:
        longest = 0
        for item in shelf:
            if (len(item) > longest):
                longest = len(item)
        
        for item in shelf:
            # If item is a long list, only print first, last element (to save screen space)
            if (type(shelf[item]) == type([]) and len(shelf[item]) > 2):
                print("{0:{arg}s} = [{1}, ..., {2}]  {3}".format(item, shelf[item][0], shelf[item][-1], type([]), arg=longest).replace('\n', ''))
            else:
                print("{0:{arg}s} = {1}  {2}".format(item, shelf[item], type(shelf[item]), arg=longest))

# For formating user messages.
def __update(tool, message):
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
                d_title = line[7:80].rstrip(); SET('d_title', d_title)

            # Get periodic box information (if any).
            if ((line[0:6]) == "CRYST1"):
                d_box   = line[7:80].rstrip(); SET('d_box', d_box)
                
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

    SET('d_residues', d_residues)

def __add_to_nameList(name):
    with shelve.open('universe') as shelf:
        # If d_nameList does not exist, add it as an empty list.
        try:
            shelf['d_nameList']
        except KeyError:
            shelf['d_nameList'] = []
    
        # Append name to d_nameList.
        temp = shelf['d_nameList']
        temp.append(name)
        shelf['d_nameList'] = temp

# Write (modified) .pdb file.
def __writepdb(name):            
    with open(name, 'w') as file:
        file.write("TITLE {0}\n".format(GET('d_title')))
        file.write("CRYST1{0}\n".format(GET('d_box')))
        file.write("MODEL {:8d}\n".format(GET('d_model')))

        atomNumber = 1
        for residue in GET('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

    __add_to_nameList(name)

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
    SET('d_fname', d_fname)
    SET('d_pdbName', d_fname[0:len(d_fname)-4])
    SET('d_model', d_model)
    SET('d_ALI', d_ALI)

    loadpdb(d_fname, d_model, d_ALI, d_chain)

    # Update d_chain from ["all"] to a list of actual chains used:
    d_chain = []                                    
    for residue in GET('d_residues'):
        d_chain.append(residue.d_chain)
    d_chain = list(set(d_chain))
    d_chain.sort()
    SET('d_chain', d_chain)

    __update("processpdb", 'filename = {0}, MODEL = {1}, ALI = {2}, chain(s) = {3}...'.format(d_fname, d_model, d_ALI, d_chain))
    
    # Some optional stuff:
    if genMissingChainsIDs:
        __update("processpdb", "generating missing chain-identifiers...")
        __genMissingChainIdentifiers()

    if resetResId:
        __update("processpdb", "resetting resid order...")
        __resetResId()

    # Write processed.pdb to file:
    __writepdb("{0}_PR1.pdb".format(GET('d_pdbName')))

    # Update index.ndx
    __generate_index()

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
        SET('d_dt', dt)
        SET('d_nsteps', nsteps)
        SET('d_nstxout', nstxout)

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
    SET('d_modelFF', d_modelFF)
    SET('d_modelWater', d_modelWater)
    SET('d_neutralTermini', d_neutralTermini)

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
            SET('d_constantpH', False)
            SET('d_restrainpH', False)

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
        os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -asp -glu -ignh -ff {2} -water {3} -ter >> builder.log 2>&1 {4}".format(GET('d_nameList')[-1], GET('d_pdbName'), d_modelFF, d_modelWater, xstr))
    else: # Still a bug here if we have multiple chains, so have to specify the termini multiple times
        os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -asp -glu -ignh -ff {2} -water {3} >> builder.log 2>&1 {4}".format(GET('d_nameList')[-1], GET('d_pdbName'), d_modelFF, d_modelWater, xstr))

    # Rebuild topology.
    __rebuild_topol()

    # To update d_residues.
    loadpdb("{0}_PR2.pdb".format(GET('d_pdbName')))

    # To update d_nameList.
    __add_to_nameList("{0}_PR2.pdb".format(GET('d_pdbName')))

    # To update index.ndx.
    __generate_index()

def protein_add_box(d_boxMargin, d_boxType='cubic'):
    __update("protein_add_box", "running gmx editconf (boxMargin = {0}, boxType = {1})...".format(d_boxMargin, d_boxType))

    os.system("gmx editconf -f {0} -o {1}_BOX.pdb -d {2} -bt {3} >> builder.log 2>&1".format(GET('d_nameList')[-1], GET('d_pdbName'), d_boxMargin, d_boxType))

    # To set d_boxMargin and d_boxType.
    SET('d_boxMargin', d_boxMargin)
    SET('d_boxType', d_boxType)

    # To update d_box.
    loadpdb("{0}_BOX.pdb".format(GET('d_pdbName')))

    # To update d_nameList.
    __add_to_nameList("{0}_BOX.pdb".format(GET('d_pdbName')))

def protein_add_buffer(bufpdbName, bufitpname, minSep=1.5):
    if not (GET('d_constantpH') and GET('d_restrainpH')):
        __update("protein_add_buffer", "skipping this step...")
        return
    
    __update("protein_add_buffer", "adding buffer molecules...")

    # Extract the required value.
    def extractMinimum():
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
        os.system("gmx insert-molecules -f {0} -o {1}_BUF.pdb -ci {2} -nmol {3} >> builder.log 2>&1".format(GET('d_nameList')[-1], GET('d_pdbName'), bufpdbName, countACID))
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

    # To update d_residues.
    loadpdb("{0}_BUF.pdb".format(GET('d_pdbName')))

    # To update d_nameList.
    __add_to_nameList("{0}_BUF.pdb".format(GET('d_pdbName')))

    # To update index.ndx.
    __generate_index()

def protein_add_water():
    __update("protein_add_water", "running gmx solvate...")
    
    os.system("gmx solvate -cp {0} -o {1}_SOL.pdb -p topol.top >> builder.log 2>&1".format(GET('d_nameList')[-1], GET('d_pdbName')))

    # To update topol.top.
    add_mol_to_topol("{0}.ff/{1}.itp".format(GET('d_modelFF'), GET('d_modelWater')), "Include water topology", "SOL", countRes('SOL'))

    # To update d_residues.
    loadpdb("{0}_SOL.pdb".format(GET('d_pdbName')))

    # To update d_nameList.
    __add_to_nameList("{0}_SOL.pdb".format(GET('d_pdbName')))

    # To update index.ndx.
    __generate_index()

def protein_add_ions():
    __update("protein_add_ions", "running gmx grompp and genion to add ions...")

    # Generate IONS.mdp (just a dummy required).
    os.system('touch IONS.mdp')

    # Add ion topology to topol.top.
    add_mol_to_topol("{0}.ff/ions.itp".format(GET('d_modelFF')), "Include ion topology")

    # Run gmx genion.
    os.system("gmx grompp -f IONS.mdp -c {0} -p topol.top -o IONS.tpr >> builder.log 2>&1".format(GET('d_nameList')[-1]))
    os.system("gmx genion -s IONS.tpr -o {0}_ION.pdb -p topol.top -pname NA -nname CL -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF".format(GET('d_pdbName')))

    # To update d_residues.
    loadpdb("{0}_ION.pdb".format(GET('d_pdbName')))

    # To update d_nameList.
    __add_to_nameList("{0}_ION.pdb".format(GET('d_pdbName')))

    # To update index.ndx.
    __generate_index()

# INDEX, JOBSCRIPT, RUN, RESET
################################################################################
def __generate_index():
    os.system("gmx make_ndx -f {0} >> builder.log 2>&1 << EOF\nq\nEOF".format(GET('d_nameList')[-1]))

def write_run(gmxPath="/usr/local/gromacs", options=""):
    __update("write_run", "gmxPath={0}, additional options= {1}".format(gmxPath, options))
    
    with open("run.sh", 'w') as file:
        file.write("#!/bin/bash\n\n")

        file.write("# Gromacs version to use:\n")
        file.write("source {0}/bin/GMXRC\n\n".format(gmxPath))

        file.write("gmx grompp -f MD.mdp -c {0} -p topol.top -n index.ndx -o MD.tpr -r {0}\n".format(GET('d_nameList')[-1]))
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

    file.write("gmx grompp -f MD.mdp -c {0} -p topol.top -n index.ndx -o MD.tpr -r {0}\n".format(GET('d_nameList')[-1]))
    
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
    generate_mdp('EM')

    __update("energy_minimize", "running gmx grompp and mdrun for energy minimization...")

    os.system("gmx grompp -f EM.mdp -c {0} -p topol.top -n index.ndx -o EM.tpr -r {0} >> builder.log 2>&1".format(GET('d_nameList')[-1]))
    os.system("gmx mdrun -s EM.tpr -o EM.trr -c {0}_EM.pdb -g EM.log -e EM.edr >> builder.log 2>&1".format(GET('d_pdbName')))

    __add_to_nameList("{0}_EM.pdb".format(GET('d_pdbName')))

def energy_tcouple():
    generate_mdp('NVT')
    
    __update("energy_tcouple", "running gmx grompp and mdrun for temperature coupling...")

    os.system("gmx grompp -f NVT.mdp -c {0} -p topol.top -n index.ndx -o NVT.tpr -r {0} >> builder.log 2>&1".format(GET('d_nameList')[-1]))
    os.system("gmx mdrun -s NVT.tpr -o NVT.trr -c {0}_NVT.pdb -g NVT.log -e NVT.edr >> builder.log 2>&1".format(GET('d_pdbName')))

    __add_to_nameList("{0}_NVT.pdb".format(GET('d_pdbName')))

def energy_pcouple():
    generate_mdp('NPT')

    __update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")

    os.system("gmx grompp -f NPT.mdp -c {0} -p topol.top -n index.ndx -o NPT.tpr -r {0} >> builder.log 2>&1".format(GET('d_nameList')[-1]))
    os.system("gmx mdrun -s NPT.tpr -o NPT.trr -c {0}_NPT.pdb -g NPT.log -e NPT.edr >> builder.log 2>&1".format(GET('d_pdbName')))    

    __add_to_nameList("{0}_NPT.pdb".format(GET('d_pdbName')))

# GENERATE CONSTANT-PH DATA
################################################################################

def generate_phdata(pH, lambdaM, nstOut, barrierE):
    # Data (hardcoded, specific for CHARMM2019)
    GLU_pKa   = 4.25
    GLU_dvdl  = [24.685, -577.05, 137.39, -172.69] # Orig Noora.
    GLU_atoms = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'] # atoms part of model
    GLU_qqA   = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    GLU_qqB   = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge
    
    ASP_pKa   = 3.65
    ASP_dvdl  = [37.822, -566.01, 117.97, -158.79] # Orig Noora.
    ASP_atoms = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'] # atoms part of model
    ASP_qqA   = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    ASP_qqB   = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge

    BUF_dvdl  = [670.1, -674.4, 83.19, -150.21] # Orig Noora.
    #           [' OW ' , ' HW1', ' HW2']
    BUF_qqA   = [-0.0656, 0.5328, 0.5328]
    BUF_qqB   = [-0.8476, 0.4238, 0.4238]

    # Skip this entire step if self.d_constantpH is false.
    if (not GET('d_constantpH')):
        __update("generate_phdata", "skipping this step...")
        return

    # Throw exception if MD.mdp does not exist.
    if (not os.path.isfile("MD.mdp")):
        raise Exception("MD.mdp does not exist! Did you generate MD.mdp before calling generate_phdata?")
    
    # Throw exception if index.ndx does not exist.
    if (not os.path.isfile("index.ndx")):
        raise Exception("index.ndx does not exist! Did you generate index.ndx before calling generate_phdata?")

    # Update user.
    __update("generate_phdata", "pH=%s, lambdaM=%s, nstOut=%s, barrierE=%s" % (pH, lambdaM, nstOut, barrierE))

    file = open('MD.mdp', 'a')

    # Formatting function.
    def addParam(name, value, comment = "NUL"):
        if (comment == "NUL"):
            file.write("{:54s} = {:13s}\n".format(name, str(value)))
        else:            
            file.write("{:54s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    file.write("\n; CONSTANT PH\n")

    # PART 1 - WRITE GENERAL PARAMETERS ########################################
    
    addParam('lambda-dynamics', 'yes')
    addParam('lambda-dynamics-simulation-ph', pH)
    addParam('lambda-dynamics-lambda-particle-mass', lambdaM)
    addParam('lambda-dynamics-update-nst', nstOut)
    addParam('lambda-dynamics-tau', 0.1) # hardcoded

    if GET('d_restrainpH'):
        addParam('lambda-dynamics-charge-constraints', 'yes')

    # Compile a list of acidic residues and their ResIDs.
    acidicResidueNameList = []; acidicResidueNumberList = []
    acidicResidueTypeList = []
    
    for residue in GET('d_residues'):
        if (residue.d_resname == 'GLU'):
            acidicResidueNameList.append('GLU')
            acidicResidueNumberList.append(residue.d_resid)
        
        if (residue.d_resname == 'ASP'):
            acidicResidueNameList.append('ASP')
            acidicResidueNumberList.append(residue.d_resid)

    if ('GLU' in acidicResidueNameList):
        acidicResidueTypeList.append('GLU')
    
    if ('ASP' in acidicResidueNameList):
        acidicResidueTypeList.append('ASP')

    if GET('d_restrainpH'):                   # If we restrain the charge 
        acidicResidueTypeList.append('BUF')   # we also have BUF.

    addParam('lambda-dynamics-number-lambda-residues', len(acidicResidueTypeList))
    
    if GET('d_restrainpH'):
        addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList) + 1)
    else:
        addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList))

    file.write('\n')

    # print(acidicResidueNameList)   # debug
    # print(acidicResidueNumberList) # debug
    # print(acidicResidueTypeList)   # debug

    # PART 2 - WRITE RESIDUE-TYPE SPECIFIC STUFF ###############################

    def writeBlock(number, name, dvdl, pKa, barrierE, qqA, qqB):

        def to_string(Input):
            string = ""
            for element in Input:
                string += str(element)
                string += " "
            return string

        addParam('lambda-dynamics-residue%s-name'              % (number), name)
        addParam('lambda-dynamics-residue%s-dvdl-coefficients' % (number), to_string(dvdl))
        addParam('lambda-dynamics-residue%s-reference-pka'     % (number), pKa)
        addParam('lambda-dynamics-residue%s-barrier'           % (number), barrierE)
        addParam('lambda-dynamics-residue%s-charges-state-A'   % (number), to_string(qqA))
        addParam('lambda-dynamics-residue%s-charges-state-B'   % (number), to_string(qqB))
        
        file.write('\n')

    for idx in range(0, len(acidicResidueTypeList)):
        if (acidicResidueTypeList[idx] == 'GLU'):
            writeBlock(idx + 1, 'GLU', GLU_dvdl, GLU_pKa, barrierE, GLU_qqA, GLU_qqB)

        if (acidicResidueTypeList[idx] == 'ASP'):
            writeBlock(idx + 1, 'ASP', ASP_dvdl, ASP_pKa, barrierE, ASP_qqA, ASP_qqB)

        if (acidicResidueTypeList[idx] == 'BUF'):
            # Multiplication is no-longer necessary because of Paul's commit on January 25th:
            # writeBlock(idx + 1, 'BUF', [i * len(acidicResidueNameList) for i in BUF_dvdl], 0, 0, BUF_qqA, BUF_qqB)
            writeBlock(idx + 1, 'BUF', BUF_dvdl, 0, 0, BUF_qqA, BUF_qqB)

    # PART 3 - WRITE INDIVIDUAL RESIDUE/LAMBDA-GROUP STUF ######################

    def writeResBlock(number, name, indexLambda, indexName):
        addParam('lambda-dynamics-atom-set%s-name'                  % (number), name)
        addParam('lambda-dynamics-atom-set%s-lambda-residues-index' % (number), indexLambda)
        addParam('lambda-dynamics-atom-set%s-index-group-name'      % (number), indexName)
        addParam('lambda-dynamics-atom-set%s-initial-lambda'        % (number), 0.5) # hardcoded
        
        if GET('d_restrainpH'):
            addParam('lambda-dynamics-atom-set%s-charge-restraint-group-index' % (number), 1)

        if (name == 'BUF'):
            addParam('lambda-dynamics-atom-set%s-buffer-residue' % (number), 'yes')
            addParam('lambda-dynamics-atom-set%s-buffer-residue-multiplier' % (number), len(acidicResidueNameList))

        file.write('\n')

    for idx in range(0, len(acidicResidueNameList)):
        writeResBlock(
                        idx + 1, 
                        acidicResidueNameList[idx],
                        acidicResidueTypeList.index(acidicResidueNameList[idx]) + 1,
                        'LAMBDA%s' % (idx + 1)
                        )

    if GET('d_restrainpH'):
        writeResBlock(
                        len(acidicResidueNameList) + 1,
                        'BUF',
                        acidicResidueTypeList.index('BUF') + 1,
                        'LAMBDA%s' % (len(acidicResidueNameList) + 1)
                        )

    file.close() # MD.mdp

    # PART 4 - APPEND THE LAMBDA INDEX GROUPS TO INDEX.NDX #####################

    file = open('index.ndx', 'a') # Append to existing index.ndx

    # Function for adding an indexList to index.ndx
    def writeTheGroup(number, indexList):
        file.write('\n[ LAMBDA{} ]\n'.format(number))
        for index in indexList:
            file.write('{} '.format(index))
        file.write('\n')

    atomCount = 1   # Keeps track of the atom number.
    grpNum    = 1   # Keeps track of the group (the LAMBDA%s).
    
    for residue in GET('d_residues'):       # loop through all residues

        indexList = []                      # clear indexList

        for atom in residue.d_atoms:        # for each residue, loop through the atoms

            if (residue.d_resname == 'GLU' and atom in GLU_atoms):
                indexList.append(atomCount)

            elif (residue.d_resname == 'ASP' and atom in ASP_atoms):
                indexList.append(atomCount)

            atomCount += 1                  # increment atomcount

        if (len(indexList) > 0):
            writeTheGroup(grpNum, indexList)
            grpNum += 1

    if GET('d_restrainpH'):

        atomCount = 1; indexList = []

        for residue in GET('d_residues'):
            for atom in residue.d_atoms:
            
                if (residue.d_resname == 'BUF'):
                    indexList.append(atomCount)
            
                atomCount += 1
        
        writeTheGroup(grpNum, indexList)

    file.close() # index.ndx

def restrain_dihedrals(resName, atomNameList, Type, phi, dphi, fc):
    __update("restrain_dihedrals", "will add restraints for {0} (all chains)...".format(resName))

    # Every chain has its own .itp file, so we loop through every file:
    for letter in GET('d_chain'):
        
        # This is to make sure we don't have multiple headers when we add multiple different restraints.
        first = False
        if not "[ dihedral_restraints ]" in open("topol_Protein_chain_{0}.itp".format(letter)).read():
            first = True

        # Append to the end of relevant .itp file:
        with open("topol_Protein_chain_{0}.itp".format(letter), 'a') as file:
            if first:
                file.write("[ dihedral_restraints ]\n")
                file.write("; ai aj ak al type phi dphi fc\n")

            # Write the atoms as an extra comment.
            file.write("; {0} {1} {2} {3}\n".format(atomNameList[0], atomNameList[1], atomNameList[2], atomNameList[3]))

            # Atomcount resets for every separate .itp file.
            count = 0
            for residue in GET('d_residues'):
                idxList = []
                for atom in residue.d_atoms:
                    # only increase atomcount when we read the relevant chain
                    if residue.d_chain == letter:
                        count += 1

                        if residue.d_resname == resName and atom in atomNameList:
                            idxList.append(count)

                # dihedrals have always four atoms so only activate if we have a length of four.                
                if len(idxList) == 4:
                    __update("restrain_dihedrals", "adding restraints for chain {0} {1}-{2}...".format(residue.d_chain, resName, residue.d_resid))
                    file.write("{:<6d} {:<6d} {:<6d} {:<6d}  {}  {}  {}  {}\n".format(idxList[0], idxList[1], idxList[2], idxList[3], Type, phi, dphi, fc))

# MAIN #########################################################################

SET('d_constantpH', True)
SET('d_restrainpH', True)

processpdb('1cvo.pdb')
write_reset()

generate_topology("charmm36-mar2019", "tip3p")

restrain_dihedrals('GLU', [' OE1', ' CD ', ' OE2', ' HE2'], 1,  0, 0, 10)
restrain_dihedrals('GLU', [' HA ', ' CA ', ' CB ', ' HB1'], 1, 60, 0, 10)

protein_add_box(d_boxMargin=1.0)
protein_add_buffer("/home/anton/GIT/phbuilder/grom/buffer.pdb", "/home/anton/GIT/phbuilder/grom/buffer.itp", 1.5)
protein_add_water()
protein_add_ions()

energy_minimize()
energy_tcouple()

generate_mdp('MD', nsteps=100000, nstxout=1000)
generate_phdata(7.0, 5.0, 1, 5.0)
write_run(gmxPath="/usr/local/gromacs_test2", options="-pme cpu")
# write_jobscript('test', 48, 1, 32, 'lindahl')

################################################################################

import matplotlib.pyplot as plt
import loaddata as load

def plotlambdacoordinates(fileName="", plotBUF=False):
    resNameList = []; resIdList = []
    for residue in GET('d_residues'):
        if (residue.d_resname == "GLU"):
            resNameList.append("GLU")
            resIdList.append(residue.d_resid)
        if (residue.d_resname == "ASP"):
            resNameList.append("ASP")
            resIdList.append(residue.d_resid)

    plt.figure()
    for idx in range(1, len(resNameList) + 1):
        t = load.Col("lambda_{0}.dat".format(idx), 1)
        x = load.Col("lambda_{0}.dat".format(idx), 2)
        
        plt.plot(t, x, label="%s-%s" % (resNameList[idx-1], resIdList[idx - 1]), linewidth=0.5)

    if (plotBUF):
        t = load.Col("lambda_{0}.dat".format(len(resNameList) + 1), 1)
        x = load.Col("lambda_{0}.dat".format(len(resNameList) + 1), 2)

        plt.plot(t, x, label="Buffer", linewidth=0.5)

    plt.xlabel("Time (ps)")
    plt.ylabel(r"$\lambda$-coordinate")

    plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
    
    plt.legend()
    plt.grid()

    if (not fileName == ""):
        plt.savefig("{0}.pdf".format(fileName))
        os.system("pdfcrop {0}.pdf {0}.pdf".format(fileName))
    else:
        plt.show()

# plotlambdacoordinates(plotBUF=True)
