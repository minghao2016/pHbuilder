import os

headRecords = [
    "HEADER", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "NUMMDL", 
    "REVDAT", "JRNL  ", "REMARK", "DBREF ", "SEQRES", "SHEET ", "SSBOND", 
    "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3",
    "SEQADV", "HET   ", "HETNAM", "HETSYN", "FORMUL", "HELIX ", "SITE  ",]

# Object for storing a single residue's data. ##################################
class Residue:
    # CONSTRUCTOR
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
    def inspect(self):
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

# Main object for reading, editing, and writing .pdb files. ####################
class PDB:
    # CONSTRUCTOR
    def __init__(self):
        self.d_model  = 1        # int       stores selected model.
        self.d_ALI    = "A"      # string    stores alternate loc. indicator.
        self.d_chainl = ["all"]  # list      stores which chain(s) to load.
        self.d_title  = ""       # string    stores TITLE line.
        
        self.d_atomLines = []   # list      stores atom lines.
        self.d_headLines = []   # list      stores header lines.
        self.d_residues  = []   # list      stores Residue objects.
    
    # Load a .pdb file.
    def loadpdb(self, fname, MODEL = 1, ALI = "A", CHAIN = ["all"]):
        self.d_model  = MODEL    # set MODEL number to load (MODEL 1 = default)
        self.d_ALI    = ALI      # set which ALI letter to load (A = default)
        self.d_chainl = CHAIN    # set which chain(s) to load (all = default)
        
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
                if (line[0:6]) in headRecords:
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
                    Residue
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
        
        print("writepdb   : info : residues=%s, MODEL=%s, ALI=%s, chain(s)=" % 
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

    # VARIOUS FUNCTIONS ########################################################

    # Print data for a specific residue.
    def inspectRes(self, resid, CHAIN = 'A'):
        for residue in self.d_residues:
            if (residue.d_resid == resid and residue.d_chain == CHAIN):
                residue.inspect()
                return
        else:
            print("inspectRes : cannot find residue with resid=%s and/or chain=%s" % (resid, CHAIN))

    # Returns total charge of protein (at neutral pH).
    def charge(self):
        charge = 0
        for residue in self.d_residues:
            if residue.d_resname in ["ARG", "LYS"]:
                charge += 1
            
            if residue.d_resname in ["ASP", "GLU"]:
                charge -= 1
        
        return charge

    # Return the title.
    def title(self):
        return self.d_title

    # Return the selected MODEL number.
    def model(self):
        return self.d_model

    # Set a (new) title.
    def setTitle(self, newTitle):
        self.d_title = newTitle

    # Print the header.
    def printHeader(self):
        for line in self.d_headLines:
            print(line, end = '')
