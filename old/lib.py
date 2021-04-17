import os
import string

class sim:
    
    def __defineGromacsCommands(self):
        # We define these because we re-use them in internal energy commands,
        # the run.sh and jobscript.sh files, and if we change something we
        # don't want to do it three times.
        self.g_EM_gromm  = "-f EM.mdp -c {0}_ION.pdb -p topol.top -n index.ndx -o EM.tpr -r {0}_ION.pdb".format(self.d_pdbName)
        self.g_NVT_gromm = "-f NVT.mdp -c {0}_EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r {0}_EM.pdb".format(self.d_pdbName)
        self.g_NPT_gromm = "-f NPT.mdp -c {0}_NVT.pdb -p topol.top -n index.ndx -o NPT.tpr -r {0}_NVT.pdb".format(self.d_pdbName)
        self.g_MD_gromm  = "-f MD.mdp -c {0}_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r {0}_NPT.pdb -maxwarn 1".format(self.d_pdbName)
        
        self.g_EM_md  = "-s EM.tpr -o EM.trr -c {0}_EM.pdb -g EM.log -e EM.edr".format(self.d_pdbName)
        self.g_NVT_md = "-s NVT.tpr -o NVT.trr -c {0}_NVT.pdb -g NVT.log -e NVT.edr".format(self.d_pdbName)
        self.g_NPT_md = "-s NPT.tpr -o NPT.trr -c {0}_NPT.pdb -g NPT.log -e NPT.edr".format(self.d_pdbName)
        self.g_MD_md  = "-v -s MD.tpr -o MD.trr -c {0}_MD.pdb -g MD.log -e MD.edr".format(self.d_pdbName)

    class __Residue: # Stores a single residue's data.
        def __init__(self, atoms, ali, resname, chain, resid, x, y, z):
            self.d_atoms   = atoms      # list      holds atom names
            self.d_ali     = ali        # list      holds alternative loc. ind.
            self.d_resname = resname    # string    holds residue name
            self.d_chain   = chain      # string    holds chain name (A, B)
            self.d_resid   = resid      # int       holds residue number
            self.d_x       = x          # list      holds x-coordiantes
            self.d_y       = y          # list      holds y-coordinates
            self.d_z       = z          # list      holds z-coordinatees

    def __init__(self):
        self.d_constantpH = True    # Turn constant-pH sim on or off.
        
        self.d_fname    = ""    # Store the name of the file that was loaded.
        self.d_model    = 0     # Store the model number.
        self.d_ALI      = ""    # Store the Alternative Location Indicator.
        self.d_chainl   = []    # Store the chain that was loaded.
        
        self.d_title    = ""    # Store the TITLE line.
        self.d_residues = []    # Store the Residue objects.

        self.d_pdbName  = ""    # Store the .pdb name (without extention).

    def processpdb(self, fname, MODEL = 1, ALI = "A", CHAIN = ["all"], resetResId = False, genMissingChains = False):
        self.loadpdb(fname, MODEL, ALI, CHAIN)    # Load the .pdb file.

        chainString = ""                          # Create userinfo message.
        if self.d_chainl[0] == "all":
            chainString = "all"
        else:
            chainString = "["
            for chain in self.d_chainl:
                chainString += chain
            chainString += "]"
                                                  # Chop off .pdb extention.
        self.d_pdbName = self.d_fname[0:len(self.d_fname)-4]
        self.__update("processpdb", "importing MODEL=%s, ALI=%s, chain(s)=%s, internal name is %s..." % (self.d_model, self.d_ALI, chainString, self.d_pdbName))

        if (genMissingChains):
            self.__update('processpdb', 'generating missing chain-identifiers...')
            self.__protein_genMissingChainIdentifiers()

        if (resetResId):
            self.__update("processpdb", "resetting resid numbering...")
            self.__protein_resetResId()

        self.__writepdb("%s_PR1.pdb" % self.d_pdbName)

        self.__defineGromacsCommands()      # Define energy-related gmx commands.

    def setconstantpH(self, value, restrain = True):
        self.d_constantpH = value
        self.d_restrainpH = restrain

    def loadpdb(self, fname, MODEL = 1, ALI = "A", CHAIN = ["all"]):
        self.d_fname    = fname         # Set internal states.
        self.d_model    = MODEL
        self.d_ALI      = ALI
        self.d_chainl   = CHAIN
        
        self.d_title    = ""            # Clear previous title.
        self.d_residues = []            # Clear previous values.

        atomLines = []
        with open(fname, 'r') as file:  # Read .pdb line-by-line.
            read = True                 # True if no specific MODEL specified.

            for line in file.readlines():
                
                # Store the .pdb TITLE in d_title.
                if ((line[0:6]) == "TITLE "):
                    self.d_title = line[6:(len(line) - 1)]

                # This is to make sure we only import the specified MODEL.
                if ((line[0:6]) == "MODEL "):
                    if ("MODEL {:8d}".format(self.d_model) in line):
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
                                atomLines.append(line)
                            # or if we want a selection of chains,
                            elif (line[21:22] in self.d_chainl):
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
                self.d_residues.append(
                    sim.__Residue
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

        atomLines.clear()                   # Clear d_atomLines to save memory.
    
    def __writepdb(self, fname):            # Write (modified) .pdb file.
        with open(fname, "w+") as file:
            file.write("TITLE %s\n" % self.d_title)          # Write title line.
            file.write("MODEL {:8d}\n".format(self.d_model)) # Write model line.

            num = 1                         # Write actual residues.
            for residue in self.d_residues:
                for idx in range(0, len(residue.d_atoms)):
                    file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', num, residue.d_atoms[idx], residue.d_ali[idx], residue.d_resname, residue.d_chain, residue.d_resid, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                    num += 1                # .pdb format string.

            file.write("TER\nENDMDL\n")     # Write EOF information.

    def __update(self, tool, message):      # For formating user messages.
        print("{:22s} : {:s}".format(tool, message))

    ############################################################################

    def protein_inspectRes(self, resid, CHAIN = 'A'):
        for residue in self.d_residues:
            if (residue.d_resid == resid and residue.d_chain == CHAIN):

                self.__update("protein_inspectRes", "resname = %s, resid = %s, chain = %s" % (residue.d_resname, residue.d_resid, residue.d_chain))

                for idx in range(0, len(residue.d_atoms)):
                    self.__update("protein_inspectRes", "{:^4s}{:1s}{:8.3f} {:8.3f} {:8.3f}".format(residue.d_atoms[idx], residue.d_ali[idx], residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))

                return

        raise Exception("Cannot find residue with resid=%s and/or chain=%s" % (resid, CHAIN))

    def protein_countRes(self, resname):    # Returns num of residues with a
        count = 0                           # specific resname.
        for residue in self.d_residues:
            if (residue.d_resname == resname):
                count += 1
        
        return count

    def __protein_resetResId(self):         # Reset the original residue numbering.
        num = 1
        for residue in self.d_residues:
            residue.d_resid = num
            num += 1

    def __protein_genMissingChainIdentifiers(self):
        # This function is based on the fact that the resid resets for every
        # separate chain. If that does not happen, this won't work.
        capIdx = 0
        caps   = string.ascii_uppercase

        for idx in range(0, len(self.d_residues)):
            try:
                current = self.d_residues[idx]
                Next    = self.d_residues[idx + 1]

                self.d_residues[idx].d_chain = caps[capIdx]

                # If the resid resets, we know we're at the end of a chain.
                if (current.d_resid + 1 != Next.d_resid):
                    capIdx += 1
            
            except IndexError:
                self.d_residues[idx].d_chain = caps[capIdx]

    ############################################################################

    def protein_add_forcefield(self, modelFF, modelWater, neutralTermini = False):
        countACID = self.protein_countRes("ASP") + self.protein_countRes("GLU")

        # USERINFO STUFF
        if (self.d_constantpH):
            self.__update("protein_add_forcefield", 'constant-pH is turned on...')
            
            if (countACID > 0):
                self.__update("protein_add_forcefield", "detected %s acidic residue(s):" % countACID)

                count = 1
                for residue in self.d_residues:
                    if (residue.d_resname in ['ASP', 'GLU']):
                        self.__update("protein_add_forcefield", "{:3s} {:<4d}".format(residue.d_resname, count))
                    count += 1
                self.__update("protein_add_forcefield", "(setting protonation state to true for all of these)")

            else:
                self.__update("protein_add_forcefield", "no acidic residues detected, turning off constant-pH...")
                self.d_constantpH = False
                self.d_restrainpH = False

        else:
            self.__update("protein_add_forcefield", 'constant-pH is turned off...')

        self.d_modelFF = modelFF        # Set internal force field variable
        self.__update("protein_add_forcefield", "using the %s force field with the %s water model..." % (modelFF, modelWater))

        # Count how many different chains we have
        listChains = []
        # atomCount  = 1
        for residue in self.d_residues:
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
        # (True) if self.d_constantpH is True. The line below sets the termini
        # to 0: NH3+, 0: COO- (charged, gmx default) if neutralTermini is False,
        # and to 1: NH2, 1: COOH (neutral) if neutralTermini is True.
        xstr = "<< EOF"
        for _ in range(0, countACID):
            xstr += "\n%d" % self.d_constantpH
        for _ in range(0, len(listChains)):
            xstr += "\n%d\n%d" % (neutralTermini, neutralTermini)
        xstr += "\nEOF"

        if (neutralTermini):
            self.__update("protein_add_forcefield", "setting termini to neutral (NH2 and COOH)...")
        else:
            self.__update("protein_add_forcefield", "setting termini to charged (NH3+ and COO-) (gmx default)...")
        
        self.__update("protein_add_forcefield", "running pdb2gmx to create topol.top...")

        # Generate topology and protonate (make neutral) all GLU and ASP:
        if (neutralTermini):
            os.system("gmx pdb2gmx -f %s_PR1.pdb -o %s_PR2.pdb -asp -glu -ignh -ff %s -water %s -ter >> builder.log 2>&1 %s" % (self.d_pdbName, self.d_pdbName, modelFF, modelWater, xstr))
        else: # Still a bug here if we have multiple chains, so have to specify the termini multiple times
            os.system("gmx pdb2gmx -f %s_PR1.pdb -o %s_PR2.pdb -asp -glu -ignh -ff %s -water %s >> builder.log 2>&1 %s" % (self.d_pdbName, self.d_pdbName, modelFF, modelWater, xstr))

        self.loadpdb("%s_PR2.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_box(self, boxSizeMargin):
        self.__update("protein_add_box", "running gmx editconf...")

        os.system("gmx editconf -f %s_PR2.pdb -o %s_BOX.pdb -c -d %s -bt cubic >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName, boxSizeMargin))

        self.loadpdb("%s_BOX.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_buffer(self, minSep):
        if (not (self.d_constantpH and self.d_restrainpH)):
            self.__update("protein_add_buffer", "skipping this step...")
            return

        self.__update("protein_add_buffer", "adding buffer molecules...")

        def extractMinimum():                   # Extract the required value
            def Float(fileName, line, col):
                for x, y in enumerate(open(fileName)):
                    if (x == line - 1):
                        return float(y.split()[col-1])
                                                    # Position of mindist in
            return Float("mindist.xvg", 25, 2)   # .xvg file.

        def clean():
            os.system("rm -f \\#mindist.xvg.*\\# \\#%s_BUF.pdb.*\\# mindist.xvg" % (self.d_pdbName))

        countACID = self.protein_countRes("ASP") + self.protein_countRes("GLU")

        attempts = 0
        while (True):           # Randomly add the buffer molecules
            os.system("gmx insert-molecules -f %s_BOX.pdb -o %s_BUF.pdb -ci buffer.pdb -nmol %s >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName, countACID))
            
                                # Determine the minimum distance between the protein and the buffers
            os.system("gmx mindist -f %s_BUF.pdb -s %s_BUF.pdb >> builder.log 2>&1 << EOF\n1\n13\nEOF" % (self.d_pdbName, self.d_pdbName))

            attempts += 1
            if (extractMinimum() >= minSep):
                break

            if (attempts > 100):
                clean()
                raise Exception("Maximum number of buffer insertion attempts exceeded (100). Try decreasing minSep or increasing boxSizeMargin.")

        self.__update("protein_add_buffer", "placing buffers took %s attempt(s) (mindist = %s)..." % (attempts, extractMinimum()))

        clean()                 # Cleanup

        self.__update("protein_add_buffer", "updating topology...")

        topList = []
        with open("topol.top", "r") as file:    # Add the buffer water's topology to our .top file:
            for line in file.readlines():       # This piece of code is kind of a hoax but it works.
                topList.append(line)

        with open("topol.top", "w+") as file:
            try:
                for idx in range(0, len(topList)):
                    file.write(topList[idx])

                    # If we see that the next line is this:
                    if (topList[idx + 1] == "; Include water topology\n"):
                        # Then insert the buffer topology before that line:
                        file.write("; Include buffer topology\n")
                        file.write("#include \"buffer.itp\"\n\n")            

            except IndexError:
                pass

            file.write("BUF\t\t\t\t\t  %s\n" % (countACID))
        topList.clear()

        self.loadpdb("%s_BUF.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_water(self):
        self.__update("protein_add_water", "running gmx solvate...")

        if (self.d_constantpH and self.d_restrainpH):
            os.system("gmx solvate -cp %s_BUF.pdb -o %s_SOL.pdb -p topol.top >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName))
        else:
            os.system("gmx solvate -cp %s_BOX.pdb -o %s_SOL.pdb -p topol.top >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName))

        self.loadpdb("%s_SOL.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_ions(self):
        self.__update("protein_add_ions", "running gmx grompp and genion to add ions...")
        
        os.system("gmx grompp -f IONS.mdp -c %s_SOL.pdb -p topol.top -o IONS.tpr >> builder.log 2>&1" % self.d_pdbName)
        os.system("gmx genion -s IONS.tpr -o %s_ION.pdb -p topol.top -pname NA -nname CL -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF" % self.d_pdbName)

        self.loadpdb("%s_ION.pdb" % self.d_pdbName) # Update internal d_residues.

    ############################################################################

    def generate_mdp(self, Type, nsteps = 25000, nstxout = 0, nstvout = 0, compress = 0, posres = False):
        self.firstLine = True

        if (Type not in ['EM', 'NVT', 'NPT', 'MD']):
            raise Exception("Unknown .mdp Type specified. Types are: EM, NVT, NPT, MD.")

        file = open("%s.mdp" % Type, 'w+')

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

        self.__update("generate_mdp", "Type=%s, nsteps=%s, nstxout=%s, nstvout=%s, compress=%s" % (Type, nsteps, nstxout, nstvout, compress))

        # POSITION RESTRAIN
        if (Type in ['EM', 'MD']):
            if (self.d_constantpH and self.d_restrainpH and posres):
                addTitle('Position restrain')
                addParam('define', '-DPOSRES -DPOSRES_BUF', 'Position restraints.')
            elif (self.d_constantpH and self.d_restrainpH):
                addTitle('Position restrain')
                addParam('define', '-DPOSRES_BUF', 'Position restraints.')

        if (Type in ['NVT', 'NPT']): # position restrain temp and press coupling
            addTitle('Position restrain')
            if (self.d_constantpH and self.d_restrainpH):
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
        if (Type == 'MD' and self.d_restrainpH):
            addParam('comm-mode', 'Linear', 'Remove center of mass translation.')
            addParam('comm-grps', 'Protein Non-Protein') # not our index but default gmx

        # OUTPUT CONTROL
        if (nstvout and compress):
            raise Exception("Cannot combine nstxout-compressed and nstvout.")
        
        addTitle("Output control")
        if (compress):
            addParam('nstxout-compressed', nstxout, 'Write frame every %.3f ps.' % (dt * nstxout))
            addParam('nstvout', nstvout, 'Write frame every %.3f ps.' % (dt * nstvout))
        else:
            addParam('nstxout', nstxout, 'Write frame every %.3f ps.' % (dt * nstxout))
            addParam('nstvout', nstvout, 'Write frame every %.3f ps.' % (dt * nstvout))

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
        
        if (self.d_modelFF[0:5].lower() == "charm"): # if we use a CHARMM force field...
            addParam('rcoulomb', 1.2, 'Berk: CHARMM was calibrated for 1.2 nm.')
            addParam('fourierspacing', 0.14, 'Berk: set this to 0.14 for CHARMM.')
        else: # Default for force fields:
            addParam('rcoulomb', 1.0, 'Coulomb cut-off (nm).')

        # VAN DER WAALS
        addTitle("Van der Waals")
        addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')

        if (self.d_modelFF[0:5].lower() == "charm"): # if we use a CHARMM force field...
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
            # This is only meaningful with integrator 'md', and we don't want
            # to use this while we have pressure coupling or doing MD.
            addTitle('Generate velocities for startup')
            addParam('gen_vel', 'yes')
            # addParam('gen_temp', 300)         # Default is also 300K.
            # addParam('gen_seed', -1)          # Default is also -1 (random).

        file.close()

    def generate_index(self):
        self.__update("generate_index", "running gmx make_ndx to create index.ndx...")
        os.system("gmx make_ndx -f {0}_ION.pdb -o index.ndx >> builder.log 2>&1 << EOF\nq\nEOF".format(self.d_pdbName))

    def generate_phdata_legacy(self, pH, lambdaM, nstOut, barrierE):
        if (not self.d_constantpH):
            self.__update("generate_phdata_legacy", "skipping this step...")
            return

        # Throw exception if MD.mdp does not exist.
        if (not os.path.isfile("MD.mdp")):
            raise Exception("MD.mdp does not exist! Did you generate MD.mdp before calling generate_phdata_legacy?")

        # Add this to end of MD.mdp (formerly done in mdp generator but now here).
        with open('MD.mdp', 'a') as file:
            file.write('\n; CONSTANT PH\n')
            file.write('lambda-dynamics      = yes           ; Enable constant pH.\n')

        self.__update("generate_phdata_legacy", "pH=%s, lambdaM=%s, nstOut=%s, barrierE=%s" % (pH, lambdaM, nstOut, barrierE))

        file = open("constant_ph_input.dat", "w+")
        
        ######### PART 1 - GENERAL INPUT SETTINGS FOR CONSTANT PH MD ###########
        def addParam(name, value): # Formatting function for parameters.
            file.write("{:21s} = {}\n".format(name, str(value)))

        addParam('ph', pH)          # Simulation pH

        countGLU  = self.protein_countRes("GLU")
        countASP  = self.protein_countRes("ASP")
        countACID = countGLU + countASP
        
        resTypeCount = 0
        if (countGLU > 0):          # If we have at least one GLU ...
            resTypeCount += 1
        if (countASP > 0):          # If we have at least one ASP ...
            resTypeCount += 1
        if (self.d_restrainpH):     # If we constrain the charge we have
            resTypeCount += 1       # at least one BUF

        addParam('nr_residues', resTypeCount) # Number of different TYPES of residues (including BUF if it exists)
        
        if (self.d_restrainpH): # Number of lambda groups (including the one for BUF if it exists)
            addParam('nr_lambdagroups', countACID + 1)
        else:  
            addParam('nr_lambdagroups', countACID)
        file.write('\n')

        addParam('m_lambda', lambdaM)               # mass of l-particles
        addParam('T_lambda', "300")                 # ref. temp. of l-particles
        addParam('tau', 0.1)                        # time constant for thermostat
        addParam('thermostat', 'v-rescale')         # 'v-rescale' or 'langevin'
        addParam('nst_lambda', nstOut)              # numSteps between output
        addParam('multistate_constraint', 'no')     # NOT RELEVANT FOR NOW
        addParam('n_multigroups', 0)                # NOT RELEVANT FOR NOW
        addParam('n_states', 1)                     # NOT RELEVANT FOR NOW
        
        if (self.d_restrainpH):
            addParam('charge_constraint', 'yes')
            addParam('N_buffers', countACID)        # Number of individual buffer molecules.
            addParam('m_buf', lambdaM)              # Mass of buffer particles.

            file.write('constrained_lambdas   = ')      
            for num in range(1, (countACID + 1) + 1):   
                file.write("%d " % num)             # Write constrained lambdas.

        else:
            addParam('charge_constraint', 'no')
            addParam('N_buffers', 0)                # Number of individual buffer molecules.
            addParam('m_buf', 0)                    # Mass of buffer particles.
            file.write('constrained_lambdas   = ')  # Write constrained lambdas.

        file.write('\n\n')

        ################ PART 2 - RESIDUE-TYPE SPECIFIC PARAMETERS #############

        def addRes1(name, n_coeffs, dvdl_coeffs, ref_pka):  # formatting function.
            addParam('residue', name)
            addParam('ref_pka', ref_pka)
            addParam('n_coeffs', n_coeffs)

            file.write("{:21s} = ".format('dvdl_coeffs'))
            for coeff in dvdl_coeffs:
                file.write('%.3f ' % (coeff))
            file.write('\n')

            file.write('\n')

        #   resName  numParams  params for ref. potential   refpKa
        if (countGLU > 0):
            addRes1('GLU', 4, [24.685, -577.05, 137.39, -172.69], 4.25) # Orig Noora.
            # addRes1('GLU', 4, [19.543, -596.473, 164.418, -178.376], 4.25) # Zondagavond.

        if (countASP > 0):
            addRes1('ASP', 4, [37.822, -566.01, 117.97, -158.79], 3.65) # Orig Noora.

        if (self.d_restrainpH): # New, but might not be necessary in newer commits.
            addRes1('BUF', 4, [i * countACID for i in [670.1, -674.4, 83.19, -150.21]], 0) # Orig Noora.
            # addRes1('BUF', 4, [i * countACID for i in [660.613, -662.534, 77.058, -134.331]], 0) # Zondagavond.

        ################## PART 3 - RESIDUE-SPECIFIC PARAMETERS ################

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

        count = 1; Lidx = 1
        for residue in self.d_residues:      # loop through all residues

            indexList = []

            if (residue.d_resname == 'ASP'):    # if we find an ASP
                addParam('name', 'ASP')                         # hardcoded
                addParam('residue_number', Lidx); Lidx += 1     # this is not related to the resid, it only has to do with the lambda.
                addParam('initial_lambda', '0.5')               # hardcoded
                addParam('barrier', barrierE)                   # parameter
                addParam('n_atoms', '5')                        # hardcoded

                for atom in residue.d_atoms:    # Add indices of relevant atoms
                    if (atom in ASP_atoms):     # of ASP to list
                        indexList.append(count)

                    count += 1

                writeIndexLine(indexList)

                for idx in range(0, len(indexList)):
                    file.write("{:<7d} {:7.3f}  {:7.3f}\n".format(
                        indexList[idx], ASP_charge1[idx], ASP_charge2[idx]))
                
                file.write('\n')

            elif (residue.d_resname == 'GLU'):
                addParam('name', 'GLU')                         # hardcoded
                addParam('residue_number', Lidx); Lidx += 1     # this is not related to the resid, it only has to do with the lambda.
                addParam('initial_lambda', '0.5')               # hardcoded
                addParam('barrier', barrierE)                   # parameter
                addParam('n_atoms', '5')                        # hardcoded

                for atom in residue.d_atoms:
                    if (atom in GLU_atoms):
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

        if (self.d_restrainpH):
            # WRITE BUFFER HEAD
            addParam('name', 'BUF')
            addParam('residue_number', Lidx)
            addParam('initial_lambda', '0.5')
            addParam('barrier', 0.0)
            addParam('n_atoms', 3 * countACID)

            # GET INDEXLIST FOR BUFFER
            indexList = []; count = 1
            
            for residue in self.d_residues:
                for atom in residue.d_atoms:

                    if (residue.d_resname == 'BUF'):
                        indexList.append(count)
                        
                    count += 1

            # WRITE INDEXLIST FOR BUFFER
            writeIndexLine(indexList)

            # WRITE CHARGES FOR BUFFER ATOMS
            for idx in range(0, len(indexList)):
                file.write("{:<7d} {:7.4f}  {:7.4f}\n".format(indexList[idx], BUF_charge1[idx % 3], BUF_charge2[idx % 3]))

        file.close()

    def generate_phdata(self, pH, lambdaM, nstOut, barrierE):
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
        BUF_qqB   = [-0.8476, 0.4238, 0.4328]

        # Skip this entire step if self.d_constantpH is false.
        if (not self.d_constantpH):
            self.__update("generate_phdata", "skipping this step...")
            return

        # Throw exception if MD.mdp does not exist.
        if (not os.path.isfile("MD.mdp")):
            raise Exception("MD.mdp does not exist! Did you generate MD.mdp before calling generate_phdata?")
        
        # Throw exception if index.ndx does not exist.
        if (not os.path.isfile("index.ndx")):
            raise Exception("index.ndx does not exist! Did you generate index.ndx before calling generate_phdata?")

        # Update user.
        self.__update("generate_phdata", "pH=%s, lambdaM=%s, nstOut=%s, barrierE=%s" % (pH, lambdaM, nstOut, barrierE))

        file = open('MD.mdp', 'a')

        # Formatting function.
        def addParam(name, value, comment = "NUL"):
            if (comment == "NUL"):
                file.write("{:54s} = {:13s}\n".format(name, str(value)))
            else:            
                file.write("{:54s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

        file.write("\n; CONSTANT PH\n")

        # PART 1 - WRITE GENERAL PARAMETERS ####################################
        
        addParam('lambda-dynamics', 'yes')
        addParam('lambda-dynamics-simulation-ph', pH)
        addParam('lambda-dynamics-lambda-particle-mass', lambdaM)
        addParam('lambda-dynamics-update-nst', nstOut)
        addParam('lambda-dynamics-tau', 0.1) # hardcoded

        if (self.d_restrainpH):
            addParam('lambda-dynamics-charge-constraints', 'yes')

        # Compile a list of acidic residues and their ResIDs.
        acidicResidueNameList = []; acidicResidueNumberList = []
        acidicResidueTypeList = []
        
        for residue in self.d_residues:
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

        if (self.d_restrainpH):                   # If we restrain the charge 
            acidicResidueTypeList.append('BUF')   # we also have BUF.

        addParam('lambda-dynamics-number-lambda-residues', len(acidicResidueTypeList))
        
        if (self.d_restrainpH):
            addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList) + 1)
        else:
            addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList))

        file.write('\n')

        # print(acidicResidueNameList)   # debug
        # print(acidicResidueNumberList) # debug
        # print(acidicResidueTypeList)   # debug

        # PART 2 - WRITE RESIDUE-TYPE SPECIFIC STUFF ###########################

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

        # PART 3 - WRITE INDIVIDUAL RESIDUE/LAMBDA-GROUP STUF ##################

        def writeResBlock(number, name, indexLambda, indexName):
            addParam('lambda-dynamics-atom-set%s-name'                  % (number), name)
            addParam('lambda-dynamics-atom-set%s-lambda-residues-index' % (number), indexLambda)
            addParam('lambda-dynamics-atom-set%s-index-group-name'      % (number), indexName)
            addParam('lambda-dynamics-atom-set%s-initial-lambda'        % (number), 0.5) # hardcoded
            
            if (self.d_restrainpH):
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

        if (self.d_restrainpH):
            writeResBlock(
                          len(acidicResidueNameList) + 1,
                          'BUF',
                          acidicResidueTypeList.index('BUF') + 1,
                          'LAMBDA%s' % (len(acidicResidueNameList) + 1)
                         )

        file.close() # MD.mdp

        # PART 4 - APPEND THE LAMBDA INDEX GROUPS TO INDEX.NDX #################

        file = open('index.ndx', 'a') # Append to existing index.ndx

        # Function for adding an indexList to index.ndx
        def writeTheGroup(number, indexList):
            file.write('\n[ LAMBDA{} ]\n'.format(number))
            for index in indexList:
                file.write('{} '.format(index))
            file.write('\n')

        atomCount = 1   # Keeps track of the atom number.
        grpNum    = 1   # Keeps track of the group (the LAMBDA%s).
        
        for residue in self.d_residues:         # loop through all residues

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

        if (self.d_restrainpH):

            atomCount = 1; indexList = []

            for residue in self.d_residues:
                for atom in residue.d_atoms:
                
                    if (residue.d_resname == 'BUF'):
                        indexList.append(atomCount)
                
                    atomCount += 1
            
            writeTheGroup(grpNum, indexList)

        file.close() # index.ndx

    ############################################################################

    def energy_minimize(self, skip = False):
        if (skip):
            self.__update("energy_minimize", "skipping energy minimization (only copying .pdb)...")
            
            os.system("cp {0}_ION.pdb {0}_EM.pdb".format(self.d_pdbName))
        else:
            self.__update("energy_minimize", "running gmx grompp and mdrun for energy minimization...")

            os.system("gmx grompp %s >> builder.log 2>&1" % self.g_EM_gromm)
            os.system("gmx mdrun  %s >> builder.log 2>&1" % self.g_EM_md)

    def energy_tcouple(self, skip = False):
        if (skip):
            self.__update("energy_tcouple", "skipping temperature coupling (only copying .pdb)...")
            self.__update("energy_tcouple", "WARNING: skipping this step means velocities will NOT be generated!")
            self.__update("energy_tcouple", "WARNING: This can have consequences for sampling of parameter space.")

            os.system("cp {0}_EM.pdb {0}_NVT.pdb".format(self.d_pdbName))
        else:
            self.__update("energy_tcouple", "running gmx grompp and mdrun for temperature coupling...")

            os.system("gmx grompp %s >> builder.log 2>&1" % self.g_NVT_gromm)
            os.system("gmx mdrun  %s >> builder.log 2>&1" % self.g_NVT_md)

    def energy_pcouple(self, skip = False):
        if (skip):
            self.__update("energy_pcouple", "skipping pressure coupling (only copying .pdb)...")
        
            os.system("cp {0}_NVT.pdb {0}_NPT.pdb".format(self.d_pdbName))
        else:
            self.__update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")
        
            os.system("gmx grompp %s >> builder.log 2>&1" % self.g_NPT_gromm)
            os.system("gmx mdrun  %s >> builder.log 2>&1" % self.g_NPT_md)

    ############################################################################

    def write_run(self, gmxPath = '/usr/local/gromacs', mode = 'default'):
        self.__update("write_run", "gmxPath={0}, mode={1}".format(gmxPath, mode))
        
        with open("run.sh", "w+") as file:
            file.write("#!/bin/bash\n\n")

            file.write("# source specified gromacs version\n")
            file.write("source %s/bin/GMXRC\n\n" % gmxPath)

            file.write("gmx grompp %s\n\n" % self.g_MD_gromm)
            
            if   (mode == 'default'):
                file.write("gmx mdrun %s\n\n" % self.g_MD_md)
            elif (mode == 'cpu'):
                file.write("gmx mdrun %s -nb cpu\n\n" % self.g_MD_md)
            elif (mode == 'gpu'):
                file.write("gmx mdrun %s -pme cpu\n\n" % self.g_MD_md)
            else:
                raise Exception("Unknown mode specified")

        os.system("chmod +x run.sh")

    def write_reset(self):
        self.__update("write_reset", "writing reset.sh...")

        with open("reset.sh", "w+") as file:
            file.write("#!/bin/bash\n\n")
            
            file.write("if [ -f \"%s_MD.pdb\" ]\nthen\n" % self.d_pdbName)
            file.write("\tread -p \"Warning: simulation has finished. Proceed? (y)\" var\n")
            file.write("else\n")
            file.write("\trm -rf \\_\\_py* charmm*\n")
            file.write("\trm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg\n")
            file.write("\trm -f \\#*\\#\n")
            file.write("\trm -f buffer.pdb %s_*.pdb\n" % self.d_pdbName)
            file.write("\trm -f run.sh reset.sh jobscript.sh\n")                   
            file.write("fi\n\n")

            file.write("if [ \"${var}\" = \"y\" ]\nthen\n")
            file.write("\trm -rf \\_\\_py* charmm*\n")
            file.write("\trm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg\n")
            file.write("\trm -f \\#*\\#\n")
            file.write("\trm -f buffer.pdb %s_*.pdb\n" % self.d_pdbName)
            file.write("\trm -f run.sh reset.sh jobscript.sh\n")            
            file.write("fi\n\n")

        os.system("chmod +x reset.sh")
    
    def write_jobscript(self, jobName, time, nodes, ntasks, queue):
        self.__update("write_jobscript", "jobname=%s, time=%s(hrs), nodes=%s, ntasks=%s, queue=%s..." % (jobName, time, nodes, ntasks, queue))

        file = open("jobscript.sh", "w+")

        def writeHead(param, value):
            file.write("#SBATCH --%s=%s\n" % (param, value))

        file.write("#!/bin/bash\n")

        writeHead("time", "%d-%.2d:00:00" % (int(time / 24), time % 24))
        writeHead("nodes", nodes)
        writeHead("ntasks", ntasks)
        writeHead("partition", queue)
        writeHead("job-name", jobName)
        writeHead("mail-user", "anton.jansen@scilifelab.se")
        writeHead("mail-type", "ALL")
        file.write('\n')

        file.write("module purge\nmodule load slurm\n")
        
        if (self.d_constantpH):
            file.write("\n# compile our custom Gromacs version on cluster backend node\n")
            file.write("mkdir build\n")
            file.write("cd build\n")
            file.write("cmake ~/gromacs_dev -DGMX_GPU=OFF -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=${PWD}/..\n")
            file.write("make -j\n")
            file.write("make install\n")
            file.write("cd ..\n")
            file.write("rm -r build\n")
            file.write("source ${PWD}/bin/GMXRC\n\n")
        else:
            file.write("module load gromacs/2020.3\n\n")

        file.write("if [ ! -f \"%s_NPT.pdb\" ]\nthen\n" % self.d_pdbName)
        
        file.write("\tgmx grompp %s\n"  % self.g_EM_gromm)
        file.write("\tgmx mdrun %s\n\n" % self.g_EM_md)

        file.write("\tgmx grompp %s\n"  % self.g_NVT_gromm)
        file.write("\tgmx mdrun %s\n\n" % self.g_NVT_md)

        file.write("\tgmx grompp %s\n"  % self.g_NPT_gromm)
        file.write("\tgmx mdrun %s\n"   % self.g_NPT_md)
        file.write("fi\n\n")

        file.write("gmx grompp %s\n"    % self.g_MD_gromm)
        
        if (self.d_constantpH):        
            file.write("gmx mdrun %s -nb cpu\n" % self.g_MD_md)
        else:
            file.write("gmx mdrun %s\n" % self.g_MD_md)

        file.close()
