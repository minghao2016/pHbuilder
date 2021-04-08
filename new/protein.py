import os

import universe
import utils
import topol

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

def process(d_fname, d_model=1, d_ALI='A', d_chain=["all"], resetResId=False):
    universe.add('d_fname', d_fname)
    universe.add('d_pdbName', d_fname[0:len(d_fname)-4])
    universe.add('d_model', d_model)
    universe.add('d_ALI', d_ALI)

    load(d_fname, d_model, d_ALI, d_chain)

    # Update d_chain from ["all"] to a list of actual chains used:
    d_chain = []                                    
    for residue in universe.get('d_residues'):
        d_chain.append(residue.d_chain)
    d_chain = list(set(d_chain))
    d_chain.sort()
    universe.add('d_chain', d_chain)

    utils.update("process", 'filename = {0}, MODEL = {1}, ALI = {2}, chain(s) = {3}...'.format(d_fname, d_model, d_ALI, d_chain))
    
    # Some optional stuff:
    if resetResId:
        utils.update("process", "resetting resid order...")
        __resetResId()

    # Write processed.pdb to file:
    __write("{0}_PR1.pdb".format(universe.get('d_pdbName')))

    # Update index.ndx
    utils.generate_index()

# Load a .pdb file into the universe.
def load(d_fname, d_model=1, d_ALI='A', d_chain=["all"]):
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
                d_title = line[7:80].rstrip(); universe.add('d_title', d_title)

            # Get periodic box information (if any).
            if ((line[0:6]) == "CRYST1"):
                d_box   = line[7:80].rstrip(); universe.add('d_box', d_box)
                
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

    universe.add('d_residues', d_residues)

# Write (modified) .pdb file.
def __write(name):            
    with open(name, 'w') as file:
        file.write("TITLE {0}\n".format(universe.get('d_title')))
        file.write("CRYST1{0}\n".format(universe.get('d_box')))
        file.write("MODEL {:8d}\n".format(universe.get('d_model')))

        atomNumber = 1
        for residue in universe.get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

    utils.add_to_nameList(name)

# Reset the original residue numbering.
def __resetResId():
    num = 1
    for residue in universe.get('d_residues'):
        residue.d_resid = num
        num += 1

# Inspect a specific residue.
def inspectRes(resid, chain='A'):
    for residue in universe.get('d_residues'):
        if (residue.d_resid == resid and residue.d_chain == chain):
            print("resname = {0}, resid = {1}, chain = {2}".format(residue.d_resname, residue.d_resid, residue.d_chain))

            for idx in range(0, len(residue.d_atoms)):
                print("{:^4s}{:1s}{:8.3f} {:8.3f} {:8.3f}".format(residue.d_atoms[idx], residue.d_ali[idx], residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))

            return

    raise Exception("Cannot find residue with resid={0} and/or chain={1}".format(resid, chain))

# Count number of residues with a specific residue name.
def countRes(resname):
    count = 0
    for residue in universe.get('d_residues'):
        if (residue.d_resname == resname):
            count += 1
    
    return count

def add_box(d_boxMargin, d_boxType='cubic'):
    utils.update("add_box", "running gmx editconf (boxMargin = {0}, boxType = {1})...".format(d_boxMargin, d_boxType))

    os.system("gmx editconf -f {0} -o {1}_BOX.pdb -d {2} -bt {3} >> builder.log 2>&1".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_boxMargin, d_boxType))

    # To set d_boxMargin and d_boxType.
    universe.add('d_boxMargin', d_boxMargin)
    universe.add('d_boxType', d_boxType)

    # To update d_box.
    load("{0}_BOX.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_BOX.pdb".format(universe.get('d_pdbName')))

def add_buffer(bufpdbName, bufitpname, minSep=1.5):
    if not (universe.get('d_constantpH') and universe.get('d_restrainpH')):
        utils.update("add_buffer", "skipping this step...")
        return
    
    utils.update("add_buffer", "adding buffer molecules...")

    # Extract the required value.
    def extractMinimum():
        def Float(fileName, line, col):
            for x, y in enumerate(open(fileName)):
                if (x == line - 1):
                    return float(y.split()[col-1])
                                                
        return Float("mindist.xvg", 25, 2) # Position of mindist in .xvg file.

    def clean():
        os.system("rm -f \\#mindist.xvg.*\\# \\#{0}_BUF.pdb.*\\# mindist.xvg".format(universe.get('d_pdbName')))

    countACID = countRes('ASP') + countRes('GLU')

    os.system("cp {0} {1}_BUF.pdb".format(universe.get('d_nameList')[-1], universe.get('d_pdbName')))
    idx = 0
    while (idx != countACID):
        # insert a buffer molecule
        os.system("gmx insert-molecules -f {0}_BUF.pdb -o {0}_BUF.pdb -ci {1} -nmol 1 >> builder.log 2>&1".format(universe.get('d_pdbName'), bufpdbName))
        os.system("rm -f \\#{0}_BUF.pdb.*\\#".format(universe.get('d_pdbName'))) # cleanup

        # get the distance of the buffer molecule to protein
        os.system("gmx mindist -f {0}_BUF.pdb -s {0}_BUF.pdb >> builder.log 2>&1 << EOF\n1\n13\nEOF".format(universe.get('d_pdbName')))
        
        if extractMinimum() >= minSep:
            idx += 1
            print(idx, end="\r")
        else:
            os.system("head -n -5 {0}_BUF.pdb > tmp.txt && mv tmp.txt {0}_BUF.pdb".format(universe.get('d_pdbName')))

    # attempts = 0
    # while (True):           # Randomly add the buffer molecules
    #     os.system("gmx insert-molecules -f {0} -o {1}_BUF.pdb -ci {2} -nmol {3} >> builder.log 2>&1".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), bufpdbName, countACID))
    #                         # Determine the minimum distance between the protein and the buffers
    #     os.system("gmx mindist -f {0}_BUF.pdb -s {0}_BUF.pdb >> builder.log 2>&1 << EOF\n1\n13\nEOF".format(universe.get('d_pdbName')))

    #     attempts += 1
    #     if (extractMinimum() >= minSep):
    #         break

    #     if (attempts > 100):
    #         clean()
    #         raise Exception("Maximum number of buffer insertion attempts exceeded (100). Try decreasing minSep or increasing boxSizeMargin.")

    # utils.update("add_buffer", "placing buffers took %s attempt(s) (mindist = %s)..." % (attempts, extractMinimum()))

    clean()

    # Add buffer topology to topol.top.
    utils.update("add_buffer", "updating topology...")
    topol.add_mol(bufitpname, "Include buffer topology", 'BUF', countACID)

    # To update d_residues.
    load("{0}_BUF.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_BUF.pdb".format(universe.get('d_pdbName')))

    # To update index.ndx.
    utils.generate_index()

def add_water():
    utils.update("add_water", "running gmx solvate...")
    
    os.system("gmx solvate -cp {0} -o {1}_SOL.pdb -p topol.top >> builder.log 2>&1".format(universe.get('d_nameList')[-1], universe.get('d_pdbName')))

    # To update topol.top.
    topol.add_mol("{0}.ff/{1}.itp".format(universe.get('d_modelFF'), universe.get('d_modelWater')), "Include water topology", "SOL", countRes('SOL'))

    # To update d_residues.
    load("{0}_SOL.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_SOL.pdb".format(universe.get('d_pdbName')))

    # To update index.ndx.
    utils.generate_index()

def add_ions():
    utils.update("add_ions", "running gmx grompp and genion to add ions...")

    # Generate IONS.mdp (just a dummy required).
    os.system('touch IONS.mdp')

    # Add ion topology to topol.top.
    topol.add_mol("{0}.ff/ions.itp".format(universe.get('d_modelFF')), "Include ion topology")

    # Run gmx genion.
    os.system("gmx grompp -f IONS.mdp -c {0} -p topol.top -o IONS.tpr >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx genion -s IONS.tpr -o {0}_ION.pdb -p topol.top -pname NA -nname CL -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF".format(universe.get('d_pdbName')))

    # To update d_residues.
    load("{0}_ION.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_ION.pdb".format(universe.get('d_pdbName')))

    # To update index.ndx.
    utils.generate_index()
