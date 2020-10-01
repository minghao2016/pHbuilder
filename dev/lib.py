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
                    if ("MODEL {:8d}".format(self.d_model)) in line:
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

            # If the resid of the next line is different, we are at the end, so
            if (self.d_atomLines[idx][22:26] != self.d_atomLines[idx + 1][22:26]):
                # put the data in a Residue object and append to d_residues:
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

    # Resets the numbering of residues. For example, we keep counting instead 
    # of starting 1 if we go from chain A to chain B.
    def resetResId(self):
        num = 1
        for residue in self.d_residues:
            residue.d_resid = num
            num += 1

    # Return the file name (with/without extensions)
    def fname(self, ext = True):
        if (ext):
            return self.d_fname
        else:
            return self.d_fname[0:len(self.d_fname)-4]

class mdpGen:
    def __init__(self, fname, Type, dt, nsteps, output, tgroups):
        self.firstLine = True
        file = open(fname, "w+")

        def addParam(name, value, comment = "NUL"):
            if (comment == "NUL"):
                file.write("{:20s} = {:13s}\n".format(name, str(value)))
            else:
                file.write("{:20s} = {:13s} ; {:13s}\n".format(name, str(value), comment))    

        def addTitle(title=""):
            if (self.firstLine):
                file.write("; %s\n" % (title.upper()))
                self.firstLine = False
            else:
                file.write("\n; %s\n" % (title.upper()))

        # PRINT UPDATE
        print("mdpGen     : writing %s..." % (fname))

        # POSITION RESTRAIN
        if (Type in ['NVT', 'NPT']): # position restrain temp and press coupling
            addTitle('Position restrain')
            addParam('define', '-DPOSRES', 'Position restrain protein.')

        # RUN CONTROL
        addTitle("Run control")

        if (Type == 'EM'): # emtol hardcored, pretty typical for normal MD.
            addParam('integrator', 'steep', 'Use steep for EM.')
            addParam('emtol', 1000, 'Stop when max force < 1000 kJ/mol/nm.')
            addParam('emstep', dt, 'Time step (ps).')

        if (Type in ['NVT', 'NPT', 'MD']):
            addParam('integrator', 'md')
            addParam('dt', dt, 'Time step (ps).')
        
        addParam('nsteps', nsteps, '%.3f ns' % ((dt * nsteps)/1000.0))

        # OUTPUT CONTROL
        addTitle("Output control")

        if (output):
            addParam('nstxout',   output, 'Write frame every %s ps.' % int(dt * output))
            addParam('nstenergy', output, 'Write frame every %s ps.' % int(dt * output))
            addParam('nstlog',    output, 'Write frame every %s ps.' % int(dt * output))

        if (not output):
            addParam('nstxout',   0, 'Only write last frame to .trr.')
            addParam('nstenergy', 0, 'Only write last frame to .edr.')
            addParam('nstlog',    0, 'Only write last frame to .log.')

        # NEIGHBOUR SEARCHING PARAMETERS
        addTitle("Neighbour searching")
        addParam('cutoff-scheme', 'Verlet', 'Related params are inferred by Gromacs.')

        # BONDED
        if (Type in ['NVT', 'NPT', 'MD']):
            addTitle("Bond parameters")
            addParam('constraints', 'h-bonds', 'Constrain h-bond vibrations.')
            addParam('constraint_algorithm', 'lincs', 'Holonomic constraints.')
            addParam('lincs_iter', 1, 'Related to accuracy of LINCS.')
            addParam('lincs_order', 4, 'Related to accuracy of LINCS.')

        # ELECTROSTATICS
        addTitle("Electrostatics")
        addParam('coulombtype', 'PME', 'Use Particle Mesh Ewald.')
        addParam('rcoulomb', 1.2, 'Berk: CHARMM was calibrated for 1.2 nm.')
        addParam('fourierspacing', 0.14, 'Berk: set this to 0.14 for CHARMM.')

        # VAN DER WAALS
        addTitle("Van der Waals")
        addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')
        addParam('rvdw', 1.2, 'Berk: CHARMM was calibrated for 1.2 nm.')
        addParam('vdw-modifier', 'force-switch', 'Berk: specific for CHARMM.')
        addParam('rvdw-switch', 1.0, 'Berk: specific for CHARMM.')

        # TEMPERATURE COUPLING
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

            if (Type == 'NPT'): # Required when restraining.
                addParam('refcoord_scaling', 'com', 'Required when -DPOSRES.')

        # PERIODIC BOUNDARY CONDITIONS
        addTitle("Periodic boundary condition")
        addParam('pbc', 'xyz', 'To keep molecule(s) in box.')

        # GENERATE VELOCITIES FOR STARTUP
        if (Type == 'NVT'):
            # This is only meaningful with integrator 'md', and we don't want
            # to use this while we have pressure coupling or doing MD.
            addTitle('Generate velocities for startup')
            addParam('gen_vel', 'yes')
            # addParam('gen_temp', 300)         # Default is also 300K.
            # addParam('gen_seed', -1)          # Default is also -1 (random).

        # ENABLE PH
        if (Type == 'MD'):
            addTitle('Constant pH')
            addParam('lambda-dynamics', 'yes', 'Enable constant pH.')

        file.close()

# Function that generates the required constant_ph_input.dat file.
################################################################################
def lambdaGen(pdbFname, pH, qqConstrain):
    protein   = PDB(pdbFname)
    countACID = protein.countRes("ASP") + protein.countRes("GLU")

    file = open("constant_ph_input.dat", "w+")

    ############ PART 1 - GENERAL INPUT SETTINGS FOR CONSTANT PH MD ############

    def addParam(name, value): # Formatting function for parameters.
        file.write("{:21s} = {:13s}\n".format(name, str(value)))

    addParam('ph', pH)                          # Simulation pH
    
    if (qqConstrain):
        addParam('nr_residues', 3)              # ASP, GLU, BUF
    
    if (not qqConstrain):
        addParam('nr_residues', 2)              # ASP, GLU

    addParam('nr_lambdagroups', countACID)      # Number of lambda groups
    file.write('\n')

    addParam('m_lambda', 10.0)                  # mass of l-particles
    addParam('T_lambda', 300)                   # ref. temp. of l-particles
    addParam('tau', 0.1)                        # time constant for thermostat
    addParam('thermostat', 'v-rescale')         # 'v-rescale' or 'langevin'
    addParam('nst_lambda', 100)                 # numSteps between output

    if (qqConstrain):
        addParam('charge_constraint', 'yes')
        addParam('N_buffers', 1)                # number of collective buffers
    else:
        addParam('charge_constraint', 'no')
        addParam('N_buffers', 0)                # number of collective buffers

    addParam('m_buf', 10.0)                     # mass of buffer particles
    addParam('multistate_constraint', 'no')     # NOT RELEVANT FOR NOW
    addParam('n_multigroups', 0)                # NOT RELEVANT FOR NOW
    file.write('\n')

    ################ PART 2 - RESIDUE-TYPE SPECIFIC PARAMETERS #################

    def addRes1(name, n_coeffs, dvdl_coeffs, ref_pka):  # formatting function.
        addParam('residue', name)
        addParam('n_coeffs', n_coeffs)

        file.write("{:21s} = ".format('dvdl_coeffs'))
        for coeff in dvdl_coeffs:
            file.write('%s ' % (coeff))
        file.write('\n')

        addParam('ref_pka', ref_pka)
        file.write('\n')

    #   resName  numParams  params for ref. potential   refpKa
    addRes1('GLU', 4, [24.685, -577.05, 137.39, -172.69], 4.25)
    addRes1('ASP', 4, [37.822, -566.01, 117.97, -158.79], 3.65)
    
    if (qqConstrain):
        addRes1('BUF', 4, [2010.3, -2023.2, 249.56, -450.63], 4.25)

    ################## PART 3 - RESIDUE-SPECIFIC PARAMETERS ####################

    ASP_atoms   = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'] # atoms part of model
    ASP_charge1 = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    ASP_charge2 = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge

    GLU_atoms   = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'] # atoms part of model
    GLU_charge1 = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    GLU_charge2 = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge

    # BUF_atoms = [' OW ' , ' HW1', ' HW2'] # atoms in buffer (water) molecule
    BUF_charge1 = [-0.0656, 0.5328, 0.5328] # protonated charge
    BUF_charge2 = [-0.8476, 0.4238, 0.4328] # deprotonated charge

    def writeIndexLine(indexList): # formatting function.
        file.write("{:21s} = ".format('index'))

        for num in indexList:
            file.write('%s ' % num)

        file.write('\n\n')

    count = 1
    for residue in protein.d_residues:      # loop through all residues

        indexList = []

        if (residue.d_resname == 'ASP'):    # if we find an ASP
            addParam('name', 'ASP')                         # hardcoded
            addParam('residue_number', residue.d_resid)     # pull from .pdb
            addParam('initial_lambda', '0.5')               # hardcoded
            addParam('barrier', 7.5)                        # parameter
            addParam('n_atoms', '5')                        # hardcoded
            
            for atom in residue.d_atoms:    # Add indices of relevant atoms
                if atom in ASP_atoms:       # of ASP to list
                    indexList.append(count)

                count += 1

            writeIndexLine(indexList)

            for idx in range(0, len(indexList)):
                file.write("{:<7d} {:7.3f}  {:7.3f}\n".format(
                    indexList[idx], ASP_charge1[idx], ASP_charge2[idx]))
            
            file.write('\n')

        elif (residue.d_resname == 'GLU'):
            addParam('name', 'GLU')                         # hardcoded
            addParam('residue_number', residue.d_resid)     # pull from .pdb
            addParam('initial_lambda', '0.5')               # hardcoded
            addParam('barrier', 7.5)                        # parameter
            addParam('n_atoms', '5')                        # hardcoded
            
            for atom in residue.d_atoms:
                if atom in GLU_atoms:
                    indexList.append(count)

                count += 1

            writeIndexLine(indexList)

            for idx in range(0, len(indexList)):
                file.write("{:<7d} {:7.3f}  {:7.3f}\n".format(
                    indexList[idx], GLU_charge1[idx], GLU_charge2[idx]))
            
            file.write('\n')

        else:
            for atom in residue.d_atoms:
                count += 1

    if (qqConstrain):
        # WRITE BUFFER HEAD
        addParam('name', 'BUF')                         # hardcoded
        addParam('residue_number', 1)                   # !!!
        addParam('initial_lambda', '0.5')               # hardcoded
        addParam('barrier', 0.0)                        # !!!
        addParam('n_atoms', 3 * countACID)              # !!!

        # GET INDEXLIST FOR BUFFER
        indexList = []; count = 1
        
        for residue in protein.d_residues:
            for atom in residue.d_atoms:

                if (residue.d_resname == 'BUF'):
                    indexList.append(count)
                    
                count += 1

        # WRITE INDEXLIST FOR BUFFER
        writeIndexLine(indexList)

        # WRITE CHARGES FOR BUFFER ATOMS
        for idx in range(0, len(indexList)):
            file.write("{:<7d} {:7.4f}  {:7.4f}\n".format(
                        indexList[idx], BUF_charge1[idx % 3], BUF_charge2[idx % 3]))

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
