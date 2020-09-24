import os

# DATA
headRecords = [
    "HEADER", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "NUMMDL", "REVDAT", 
    "JRNL  ", "REMARK", "DBREF ", "SEQRES", "SHEET ", "SSBOND", "CRYST1",
    "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3"]

class Residue:
    # CONSTRUCTOR
    def __init__(self, atoms, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms      # list
        self.d_resname = resname    # string
        self.d_chain   = chain      # string
        self.d_resid   = resid      # int
        self.d_x       = x          # list
        self.d_y       = y          # list
        self.d_z       = z          # list

    # DEBUG / Print Residue data.
    def inspect(self):
        print ("resname = %s, resid = %s, chain = %s" % 
              (self.d_resname, self.d_resid, self.d_chain))
            
        for idx in range(0, len(self.d_atoms)):
            print ("%s \t %.3f \t %.3f \t %.3f" %
                  (self.d_atoms[idx], self.d_x[idx], self.d_y[idx], 
                   self.d_z[idx]))

class PDB:
    # CONSTRUCTOR
    def __init__(self):
        self.d_model = 1        # int       stores selected model.
        self.d_read  = False    # bool      stores bool for reading atoms.
        self.d_title = ""       # string    stores TITLE line.
        
        self.d_atomLines = []   # list      stores atom lines.
        self.d_headLines = []   # list      stores header lines.
        self.d_residues  = []   # list      stores Residue objects.
    
    # Load a .pdb file.
    def load(self, fname, MODEL = 1):
        
        # Possibly, we have multiple MODELs in the .pdb file.
        # Default is we only load MODEL 1 unless otherwise specified.
        self.d_model = MODEL
        print("Importing MODEL %s from %s..." % (self.d_model, fname))

        # Read .pdb line by line
        with open(fname, 'r') as file:
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
                        self.d_read = True
                    else:
                        self.d_read = False

                # Import atoms for specified MODEL.
                if (line[0:6]) == "ATOM  " and self.d_read:
                    self.d_atomLines.append(line)

        # Add one line of padding to prevent IndexError
        self.d_atomLines.append("000000000000000000000000000000000000000000000")

        # Loop through lines and create list of Residue objects.
        atoms = []; xCoord = []; yCoord = []; zCoord = []
        
        for idx in range(0, len(self.d_atomLines) - 1):
            atoms.append(self.d_atomLines[idx][12:16])
            xCoord.append(float(self.d_atomLines[idx][30:38]))
            yCoord.append(float(self.d_atomLines[idx][38:46]))
            zCoord.append(float(self.d_atomLines[idx][46:54]))

            if (self.d_atomLines[idx][22:26] != self.d_atomLines[idx + 1][22:26]):

                self.d_residues.append(
                    Residue
                    (
                        atoms, 
                        self.d_atomLines[idx][17:20], 
                        self.d_atomLines[idx][21:22], 
                        int(self.d_atomLines[idx][22:26]), 
                        xCoord, 
                        yCoord,
                        zCoord
                    ))

                # Reset.
                atoms = []; xCoord = []; yCoord = []; zCoord = []

        # Clear d_atomLines to save memory.
        self.d_atomLines.clear()

    # Write (modified) .pdb file.
    def write(self, fname, HEADER = False):
        
        # Give user some information about what we are doing
        if (HEADER):
            print ("Writing data (MODEL %s) to %s (including header)..." % 
                  (self.d_model, fname))
        else:
            print ("Writing data (MODEL %s) to %s (excluding header)..." %
                  (self.d_model, fname))

        file = open(fname, "w+")
        
        # Write title.
        file.write("TITLE     %s\n" % self.d_title)

        # Write the header lines if HEADER is True.
        if (HEADER):
            for line in self.d_headLines:
                file.write(line)

        # Write the model line.
        file.write("MODEL        %s\n" % self.d_model)

        # Write actual residues.
        num = 1
        for residue in self.d_residues:
            for idx in range(0, len(residue.d_atoms)):
                file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', num, residue.d_atoms[idx], ' ', residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))

                num += 1

        # Write EOF information.
        file.write("TER\nENDMDL\n")

        # Close flile.
        file.close()

    # Print data for a specific residue.
    def inspectRes(self, resid):
        resid -= 1
        self.d_residues[resid].inspect()

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

    # Return the MODEL number.
    def model(self):
        return self.d_model

    # Set a (new) title.
    def setTitle(self, newTitle):
        self.d_title = newTitle

    # Print the header.
    def printHeader(self):
        for line in self.d_headLines:
            print(line, end = '')
