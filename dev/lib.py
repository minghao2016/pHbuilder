import os

# For doing all kinds of stuff related to .pdb files.
class PDB:
    # Data.
    headRecords = [ 
        "HEADER", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "NUMMDL", 
        "REVDAT", "JRNL  ", "REMARK", "DBREF ", "SEQRES", "SHEET ", "SSBOND", 
        "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3",
        "SEQADV", "HET   ", "HETNAM", "HETSYN", "FORMUL", "HELIX ", "SITE  ",]

     # Stores a single residue's data.
    class Residue:
        def __init__(self, atoms, ali, resname, chain, resid, x, y, z):
            self.d_atoms   = atoms      # list
            self.d_ali     = ali        # list
            self.d_resname = resname    # string
            self.d_chain   = chain      # string
            self.d_resid   = resid      # int
            self.d_x       = x          # list
            self.d_y       = y          # list
            self.d_z       = z          # list

        # Print Residue data.
        def __inspect(self):
            print ("inspectRes : resname = %s, resid = %s, chain = %s" % 
                (self.d_resname, self.d_resid, self.d_chain))
                
            for idx in range(0, len(self.d_atoms)):
                print ("           :{:^4s}{:1s}{:8.3f} {:8.3f} {:8.3f}".format(
                        self.d_atoms[idx], 
                        self.d_ali[idx], 
                        self.d_x[idx], 
                        self.d_y[idx], 
                        self.d_z[idx]
                    ))

    def __init__(self, fname, MODEL = 1, ALI = "A", CHAIN = ["all"]):
        self.d_fname  = fname    # Store the name of the file that was loaded
        self.d_model  = MODEL    # set MODEL number to load (MODEL 1 = default)
        self.d_ALI    = ALI      # set which ALI letter to load (A = default)
        self.d_chainl = CHAIN    # set which chain(s) to load (all = default)

        self.d_title     = ""    # string    stores TITLE line.
        self.d_atomLines = []    # list      stores atom lines.
        self.d_headLines = []    # list      stores header lines.
        self.d_residues  = []    # list      stores Residue objects.
        
        # PRINT USER INFO ######################################################
        print("loadpdb    : importing MODEL=%s, ALI=%s, chain(s)=" % (
            self.d_model, self.d_ALI), end = '')
        
        if self.d_chainl[0] == "all":
            print("all", end = '')
        else: 
            print("[", end = '')
            for chain in self.d_chainl: 
                print("%s" % (chain), end = '')
            print("]", end = '')

        print(" from \"%s\"..." % (fname))

        # Read .pdb line by line ###############################################
        with open(fname, 'r') as file:
            read = True         # Needs to be True if no MODEL is specified.

            for line in file.readlines():
                
                # Store the .pdb TITLE in d_title.
                if (line[0:6]) == "TITLE ":
                    self.d_title = line[10:(len(line) - 1)]

                # Store all header lines in d_headLines.
                if (line[0:6]) in PDB.headRecords:
                    self.d_headLines.append(line)

                # This is to make sure we only import the specified MODEL.
                if (line[0:6]) == "MODEL ":
                    if str(self.d_model) in line:
                        read = True
                    else:
                        read = False
                
                # if our line specifies an ATOM,
                if (line[0:6] == "ATOM  "):                 
                    # and we are currently reading the correct MODEL,
                    if (read == True):
                        # and our line contains the correct specified alternate 
                        # location specifier (if any at all),
                        if (line[16:17] in [ALI, " "]):
                            # and we want all chains,
                            if (self.d_chainl == ["all"]):
                                # then load the atom.
                                self.d_atomLines.append(line)
                            # or if we want a selection of chains,
                            elif (line[21:22] in self.d_chainl):
                                # load that selection.
                                self.d_atomLines.append(line)

        # Add one line of padding to prevent IndexError
        self.d_atomLines.append("000000000000000000000000000000000000000000000")

        # Loop through lines and create list of Residue objects.
        atoms = []; ali = []; xCoord = []; yCoord = []; zCoord = []
        
        for idx in range(0, len(self.d_atomLines) - 1):
            atoms.append(self.d_atomLines[idx][12:16])
            ali.append(self.d_atomLines[idx][16:17])
            xCoord.append(float(self.d_atomLines[idx][30:38]))
            yCoord.append(float(self.d_atomLines[idx][38:46]))
            zCoord.append(float(self.d_atomLines[idx][46:54]))

            if (self.d_atomLines[idx][22:26] != self.d_atomLines[idx + 1][22:26]):

                self.d_residues.append(
                    PDB.Residue
                    (
                        atoms, 
                        ali,
                        self.d_atomLines[idx][17:20], 
                        self.d_atomLines[idx][21:22], 
                        int(self.d_atomLines[idx][22:26]), 
                        xCoord, 
                        yCoord,
                        zCoord
                    ))

                # Reset.
                atoms = []; ali = []; xCoord = []; yCoord = []; zCoord = []

        # Clear d_atomLines to save memory.
        self.d_atomLines.clear()

    # Write (modified) .pdb file.
    def writepdb(self, fname, HEADER = False):
        
        # PRINT USER INFO ######################################################
        print("writepdb   : writing \"%s\" to \"%s\"..." % (self.d_title, fname))

        print("           : info : residues=%s, MODEL=%s, ALI=%s, chain(s)=" % 
             (len(self.d_residues), self.d_model, self.d_ALI), end='')
        
        if self.d_chainl[0] == "all":
            print("all", end = '')
        else: 
            print("[", end = '')
            for chain in self.d_chainl: 
                print("%s" % (chain), end = '')
            print("]", end = '')        
        
        print(", header=%s" % (HEADER))

        # WRITE ################################################################
        with open(fname, "w+") as file:
            # write title
            file.write("TITLE     %s\n" % self.d_title)     

            # write header if HEADER
            if (HEADER):                             
                for line in self.d_headLines:
                    file.write(line)

            # write model line
            file.write("MODEL        %s\n" % self.d_model)

            # Write actual residues.
            num = 1
            for residue in self.d_residues:
                for idx in range(0, len(residue.d_atoms)):
                    file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', num, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                    num += 1

            # Write EOF information.
            file.write("TER\nENDMDL\n")

    # WRITE INDEX FILE #########################################################
    
    def writendx(self, fname, name, group):
        # Warn user if the energy group already exists
        if (os.path.isfile(fname)):
            print("writendx   : Detected existing file %s, will append [ %s ]..." % (fname, name))

            with open("index.ndx", "r") as file:
                if name in file.read():
                    print("           : Warning : [ %s ] in %s already exists. aborting..." % (name, fname))
                    return
        else:
            print("writendx   : No existing file %s was found. Will create..." % (fname))

        with open("index.ndx", "a+") as file:
            file.write("[ %s ]\n" % (name))

            count = 0; atom = 1; xxx = False
            for residue in self.d_residues:
                for _ in range(0, len(residue.d_atoms)):
                
                    if residue.d_resname in group:
                        file.write("{:<6d} ".format(atom))
                        count += 1
                        xxx = True

                        if (count == 11): # Keep rows within 80 chars.
                            file.write("\n")
                            count = 0

                    atom += 1

            file.write("\n\n")
        
        if (xxx):
            print("           : Wrote group [ %s ] from %s to %s" % (name, self.d_fname, fname))
        else:
            print("           : Warning : no atoms beloning to [ %s ] were found" % (name))

    # VARIOUS FUNCTIONS ########################################################

    # Print data for a specific residue.
    def inspectRes(self, resid, CHAIN = 'A'):
        for residue in self.d_residues:
            if (residue.d_resid == resid and residue.d_chain == CHAIN):
                residue.__inspect()
                return
        else:
            print("inspectRes : Cannot find residue with resid=%s and/or chain=%s" % (resid, CHAIN))

    # Returns the number of residues of a specific type in the protein.
    def countRes(self, resname):
        count = 0
        
        for residue in self.d_residues:
            if residue.d_resname == resname:
                count += 1
        
        return count

    # Return the file name (with/without extensions)
    def fname(self, ext = True):
        if (ext):
            return self.d_fname
        else:
            return self.d_fname[0:len(self.d_fname)-4]

def addParam(fname, name, value, comment = "NUL"):
    with open(fname, "a+") as file:
        if (comment == "NUL"):
            file.write("{:16s} = {:13s}\n".format(name, value))
        else:
            file.write("{:16s} = {:13s} ; {:13s}\n".format(name, value, comment))

class mdpGen:
    def __init__(self, fname):
        self.firstLine = True
        file  = open(fname, "w+")

        def addParam(name, value, comment = "NUL"):
            if (comment == "NUL"):
                file.write("{:16s} = {:13s}\n".format(name, str(value)))
            else:
                file.write("{:16s} = {:13s} ; {:13s}\n".format(name, str(value), comment))    

        def addTitle(title=""):
            if (self.firstLine):
                file.write("; %s\n" % (title.upper()))
                self.firstLine = False
            else:
                file.write("\n; %s\n" % (title.upper()))

        # Run control parameters
        addTitle("Run control")

        if (fname == "EM.mdp"):
            addParam('integrator', 'steep', 'Use steep for EM')
            addParam('emtol', 1000)
            addParam('emstep', 0.0005)
            
        else:
            addParam('integrator', 'md')
            addParam('dt', 0.002, 'time step (ps)')
            addParam('nsteps', 500000, 'dt * nsteps = 1 ns')
        
        # addParam('comm-mode', 'Linear')       # find out of this is necessary
        # addParam('nstcomm', 10)               # find out of this is necessary

        # Output control
        addTitle("Output control")
        addParam('nstxout', 5000, 'output one frame every ... ns')

        # Bond parameters       # is this whole section even necessary?
        addTitle("Bonds")
        addParam('constraint_algorithm', 'lincs', 'holotomic constraints')

        # Neighbour searching parameters
        addTitle("Neighbour searching")
        addParam('cutoff-scheme', 'Verlet')
        addParam('nstlist', 10)
        # addParam('ns_type', 'grid', 'neighbour searching algorithm (simple or grid)')
        # Apparently ns_type is obsolete
        addParam('rlist', 1.0, 'nblist cut-off ')

        # Electrostatics parameters
        addTitle("Electrostatics")
        addParam('coulombtype', 'PME')
        addParam('rcoulomb', 1.0, 'coulomb cutoff radius (nm)')
        addParam('epsilon_r', 80, 'relative dielectric constant for the medium')

        # Van der Waals
        addTitle("Van der Waals")
        addParam('vdwtype', 'cut-off')
        addParam('rvdw', 1.0)

        if (not fname == "EM.mdp"):  # No temperature or pressure coupling when
            # Temperature coupling      we are running EM
            addTitle("Temperature coupling")
            addParam('tcoupl', 'v-rescale')
            addParam('tc-grps', 'PROTEIN WATER_IONS', 'groups to couple seperately')
            addParam('tau-t', '0.1 0.1', 'time constant')
            addParam('ref-t', '300 3000', 'reference temperature (K)')

        # Pressure coupling
            addTitle('Pressure coupling')
            
            addParam('pcoupl', 'parrinello-rahman')
            addParam('pcoupltype', 'isotropic', 'uniform scaling of box')
            addParam('tau_p', '2.0', 'time constant (ps)')
            addParam('ref_p', 1.0, 'reference pressure (bar)')
            addParam('compressibility', 4.5e-5, 'isothermal compressbility of water')

        # Lambda parameters
        addTitle("Lambda/constant-pH")

        # Periodic boundary conditions
        addTitle("Periodic boundary condition")
        addParam('pbc', 'xyz', 'to keep molecule(s) in box')

        # Generate velocities for startup
        addTitle('Generate velocities for startup')
        addParam('gen_vel', 'yes')
        addParam('gen_temp', 300)
        addParam('gen_seed', -1)

        file.close()

# LOOSE FUNCTIONS ##############################################################

def backupFile(fname):
    if os.path.isfile(fname):
        num = 1
        while (True):
            if os.path.isfile("#%s.%s#" % (fname, num)):
                num += 1
                continue

            os.system("mv %s '#%s.%s#'" % (fname, fname, num))
            break
