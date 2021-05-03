import os
import universe, utils, topol

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

def process(fname, d_model=1, d_ALI='A', d_chain=["all"], resetResId=False):
    basename = os.path.basename(fname)
    universe.add('d_pdbName', basename[0:len(basename)-4])
    universe.add('d_model', d_model)
    universe.add('d_ALI', d_ALI)

    load(fname, d_model, d_ALI, d_chain)

    # Update d_chain from ["all"] to a list of actual chains used:
    d_chain = []                                    
    for residue in universe.get('d_residues'):
        d_chain.append(residue.d_chain)
    d_chain = list(set(d_chain))
    d_chain.sort()
    universe.add('d_chain', d_chain)

    utils.update("process", 'file = {0}, MODEL# = {1}, ALI = {2}, chain(s) = {3}...'.format(fname, d_model, d_ALI, d_chain))

    # Write processed.pdb to file:
    write("{0}_PR1.pdb".format(universe.get('d_pdbName')))

# Load a .pdb file into the universe.
def load(fname, d_model=1, d_ALI='A', d_chain=["all"]):
    d_residues = []
    atomLines  = []
    with open(fname) as file:       # Read .pdb line-by-line.
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
def write(name):            
    with open(name, 'w') as file:
        file.write("TITLE {0}\n".format(universe.get('d_title')))
        
        if universe.has('d_box'):
            file.write("CRYST1{0}\n".format(universe.get('d_box')))
        
        file.write("MODEL {:8d}\n".format(universe.get('d_model')))

        atomNumber = 1
        for residue in universe.get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

    utils.add_to_nameList(name)

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

def add_buffer(ph_bufpdbName, ph_bufitpName, ph_bufMargin=2.0, ph_bufnmol=-1, attempts=10000):
    if not (universe.get('ph_constantpH') and universe.get('ph_restrainpH')):
        utils.update("add_buffer", "either ph_constantpH or ph_restrainpH is False --> skipping...")
        return

    # If user doesn't specified the amount, use #BUF = #ACID.
    if (ph_bufnmol == -1):
        ph_bufnmol = countRes('ASP') + countRes('GLU')

    utils.update("add_buffer", "will attempt to add {0} buffer molecules...".format(ph_bufnmol))

    # RUN GROMACS INSERT-MOLECULES COMMAND
    os.system("touch vdwradii.dat") # we need this dummy file for this to work.

    os.system("gmx insert-molecules -f {0} -o {1}_BUF.pdb -ci {2} -nmol {3} -scale 1.0 -radius {4} -try {5} >> builder.log 2>&1".format(
        universe.get('d_nameList')[-1],
        universe.get('d_pdbName'),
        ph_bufpdbName,
        ph_bufnmol,
        0.5 * ph_bufMargin,
        int(attempts / ph_bufnmol)))

    os.remove("vdwradii.dat") # clean dummy file.

    # To update d_residues.
    load("{0}_BUF.pdb".format(universe.get('d_pdbName')))    

    # Give user a warning if there wasn't enough space.
    actual = countRes('BUF')
    if actual < ph_bufnmol:
        utils.update("add_buffer", "warning: only {0}/{1} requested buffer molecules inserted after {2} attempts,".format(actual, ph_bufnmol, attempts))
        utils.update("add_buffer", "warning: try decreasing ph_bufMargin (={0}nm) or increasing d_boxMargin (={1}nm)...".format(ph_bufMargin, universe.get('d_boxMargin')))
    else:
        utils.update("add_buffer", "succesfully added {0} buffer molecules...".format(actual))

    # To add buffer topology to topol.top.
    utils.update("add_buffer", "updating topology...")
    os.system("cp {} .".format(ph_bufitpName))
    topol.add_mol(os.path.basename(ph_bufitpName), "Include buffer topology", 'BUF', actual)

    # Set some parameters in the universe.
    universe.add('ph_bufpdbName', ph_bufpdbName)
    universe.add('ph_bufitpName', ph_bufitpName)
    universe.add('ph_bufMargin', ph_bufMargin)
    universe.add('ph_bufnmol', actual)

    # To update d_nameList.
    utils.add_to_nameList("{0}_BUF.pdb".format(universe.get('d_pdbName')))

def add_water():
    utils.update("add_water", "running gmx solvate...")
    
    os.system("gmx solvate -cp {0} -o {1}_SOL.pdb >> builder.log 2>&1".format(universe.get('d_nameList')[-1], universe.get('d_pdbName')))

    # To update d_residues.
    load("{0}_SOL.pdb".format(universe.get('d_pdbName')))

    # To update topol.top.
    topol.add_mol("{0}.ff/{1}.itp".format(universe.get('d_modelFF'), universe.get('d_modelWater')), "Include water topology", "SOL", countRes('SOL'))

    # To update d_nameList.
    utils.add_to_nameList("{0}_SOL.pdb".format(universe.get('d_pdbName')))

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
