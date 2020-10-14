import os

class sim:
    
    def __defineGromacsCommands(self):
        # We define these because we re-use them in internal energy commands,
        # the run.sh and jobscript.sh files, and if we change something we
        # don't want to do it three times.
        self.g_EM_gromm  = "-f EM.mdp -c {0}_ION.pdb -p topol.top -n index.ndx -o EM.tpr -r {0}_ION.pdb".format(self.d_pdbName)
        self.g_NVT_gromm = "-f NVT.mdp -c {0}_EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r {0}_EM.pdb".format(self.d_pdbName)
        self.g_NPT_gromm = "-f NPT.mdp -c {0}_NVT.pdb -p topol.top -n index.ndx -o NPT.tpr -r {0}_NVT.pdb".format(self.d_pdbName)
        self.g_MD_gromm  = "-f MD.mdp -c {0}_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r {0}_NPT.pdb".format(self.d_pdbName)
        
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

    def processpdb(self, fname, MODEL = 1, ALI = "A", CHAIN = ["all"]):
        self.__loadpdb(fname, MODEL, ALI, CHAIN)    # Load the .pdb file.

        chainString = ""                            # Create userinfo message.
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
        
        self.__writepdb("%s_PR1.pdb" % self.d_pdbName)

        self.__defineGromacsCommands()      # Define energy-related gmx commands.

    def setconstantpH(self, value):
        self.d_constantpH = value

    def __loadpdb(self, fname, MODEL = 1, ALI = "A", CHAIN = ["all"]):
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
                if (line[0:6]) == "TITLE ":
                    self.d_title = line[10:(len(line) - 1)]

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

################################################################################

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
            if residue.d_resname == resname:
                count += 1
        
        return count

    def protein_resetResId(self):           # Reset the original residue numbering.
        num = 1
        for residue in self.d_residues:
            residue.d_resid = num
            num += 1

################################################################################

    def protein_add_forcefield(self, modelFF, modelWater):
        # User update
        countACID = self.protein_countRes("ASP") + self.protein_countRes("GLU")
        self.__update("protein_add_forcefield", "detected %s acidic residues:" % countACID)

        count = 1
        for residue in self.d_residues:
            if residue.d_resname in ['ASP', 'GLU']:
                self.__update("protein_add_forcefield", "{:3s} {:<4d}".format(residue.d_resname, count))
            count += 1

        self.d_modelFF = modelFF        # Set internal force field variable
        self.__update("protein_add_forcefield", "using the %s force field with the %s water model..." % (modelFF, modelWater))

        xstr = "<< EOF"                 # Create EOF string required for pdb2gmx 
        for _ in range(0, countACID):   # to set the protonation state of ASP 
            xstr += "\n%d" % self.d_constantpH
        xstr += "\nEOF"                 # and GLU to true (specify 1 for user input option.

        self.__update("protein_add_forcefield", "running pdb2gmx to create topol.top...")

        # Generate topology and protonate (make neutral) all GLU and ASP:
        os.system("gmx pdb2gmx -f %s_PR1.pdb -o %s_PR2.pdb -asp -glu -ignh -ff %s -water %s >> builder.log 2>&1 %s" % (self.d_pdbName, self.d_pdbName, modelFF, modelWater, xstr))

        self.__loadpdb("%s_PR2.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_box(self, boxSizeMargin):
        self.__update("protein_add_box", "running gmx editconf...")

        os.system("gmx editconf -f %s_PR2.pdb -o %s_BOX.pdb -c -d %s -bt cubic >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName, boxSizeMargin))

        self.__loadpdb("%s_BOX.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_buffer(self, minSep):
        if (not self.d_constantpH):
            self.__update("protein_add_buffer", "skipping this step...")
            return

        self.__update("protein_add_buffer", "adding buffer molecules...")

        def extractMinimum():                   # Extract the required value
            def getfloat(file, line, column):   # from our mindist.xvg
	            x = open(file,'r')
	            for y, z in enumerate(x):
		            if y == line-1:
			            number = z.split()[column-1]
	            x.close()
	            return float(number)
                                                    # Position of mindist in
            return getfloat("mindist.xvg", 25, 2)   # .xvg file.

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
                    if topList[idx + 1] == "; Include water topology\n":
                        # Then insert the buffer topology before that line:
                        file.write("; Include buffer topology\n")
                        file.write("#include \"buffer.itp\"\n\n")            

            except IndexError:
                pass

            file.write("BUF\t\t\t\t\t  %s\n" % (countACID))
        topList.clear()

        self.__loadpdb("%s_BUF.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_water(self):
        self.__update("protein_add_water", "running gmx solvate...")

        if (self.d_constantpH):
            os.system("gmx solvate -cp %s_BUF.pdb -o %s_SOL.pdb -p topol.top >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName))
        else:
            os.system("gmx solvate -cp %s_BOX.pdb -o %s_SOL.pdb -p topol.top >> builder.log 2>&1" % (self.d_pdbName, self.d_pdbName))

        self.__loadpdb("%s_SOL.pdb" % self.d_pdbName) # Update internal d_residues.

    def protein_add_ions(self):
        self.__update("protein_add_ions", "running gmx grompp and genion to add ions...")
        
        os.system("gmx grompp -f IONS.mdp -c %s_SOL.pdb -p topol.top -o IONS.tpr >> builder.log 2>&1" % self.d_pdbName)
        os.system("gmx genion -s IONS.tpr -o %s_ION.pdb -p topol.top -pname NA -nname CL -neutral >> builder.log 2>&1 << EOF\nSOL\nEOF" % self.d_pdbName)

        self.__loadpdb("%s_ION.pdb" % self.d_pdbName) # Update internal d_residues.

################################################################################

    def generate_mdp(self, Type, dt, nsteps, output, tgroups):
        self.firstLine = True

        if Type not in ['EM', 'NVT', 'NPT', 'MD']:
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

        self.__update("generate_mdp", "Type=%s, dt=%s, nsteps=%s, output=%s" % (Type, dt, nsteps, output))

        # POSITION RESTRAIN
        if (Type in ['EM', 'MD'] and self.d_constantpH):
            addTitle('Position restrain')
            addParam('define', '-DPOSRES_BUF', 'Position restraints.')

        if (Type in ['NVT', 'NPT']): # position restrain temp and press coupling
            if (self.d_constantpH):
                addParam('define', '-DPOSRES -DPOSRES_BUF', 'Position restraints.')
            else:
                addTitle('Position restrain')
                addParam('define', '-DPOSRES', 'Position restraints.')

        # RUN CONTROL
        addTitle("Run control")

        if (Type == 'EM'): # emtol hardcored, pretty typical for normal MD.
            addParam('integrator', 'steep', 'Use steep for EM.')
            addParam('emtol', 1000, 'Stop when max force < 1000 kJ/mol/nm.')
            addParam('emstep', dt, 'Time step (ps).')

        if (Type in ['NVT', 'NPT', 'MD']):
            addParam('integrator', 'md')
            addParam('dt', dt, 'Time step (ps).')
        
        addParam('nsteps', nsteps, '%.1f ns' % ((dt * nsteps)/1000.0))

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

        # ENABLE PH
        if (Type == 'MD' and self.d_constantpH):
            addTitle('Constant pH')
            addParam('lambda-dynamics', 'yes', 'Enable constant pH.')

        file.close()

    def generate_index(self, name, group):
        # Warn user if the energy group already exists
        if (os.path.isfile("index.ndx")):
            self.__update("generate_index", "detected existing index.ndx, appending [ %s ]..." % name)

            with open("index.ndx", "r") as file:
                if name in file.read():
                    self.__update("generate_index", "warning : [ %s ] in index.ndx already exists. Skipping..." % name)
                    return
        else:
            self.__update("generate_index", "no existing index.ndx found. Will create...")
            self.__update("generate_index", "writing [ %s ]..." % name)

        with open("index.ndx", "a+") as file:
            file.write("[ %s ]\n" % (name))

            count = 0; atom = 1; foundAtLeastOneAtom = False
            for residue in self.d_residues:
                for _ in range(0, len(residue.d_atoms)):
                
                    if residue.d_resname in group:
                        file.write("{:<6d} ".format(atom))
                        count += 1
                        foundAtLeastOneAtom = True

                        if (count == 11): # Keep rows within 80 chars.
                            file.write("\n")
                            count = 0

                    atom += 1

            file.write("\n\n")
        
        if (not foundAtLeastOneAtom):
            self.__update("generate_index", "no atoms belonging to [ %s ] were found..." % name)

    def generate_phdata(self, pH):
        if (not self.d_constantpH):
            self.__update("generate_phdata", "skipping this step...")
            return

        self.__update("generate_phdata", "pH=%s" % pH)
        
        countACID = self.protein_countRes("ASP") + self.protein_countRes("GLU")

        file = open("constant_ph_input.dat", "w+")

        ######### PART 1 - GENERAL INPUT SETTINGS FOR CONSTANT PH MD ###########

        def addParam(name, value): # Formatting function for parameters.
            file.write("{:21s} = {:13s}\n".format(name, str(value)))

        addParam('ph', pH)                          # Simulation pH
        addParam('nr_residues', 3)                  # ASP, GLU, BUF
        addParam('nr_lambdagroups', countACID)      # Number of lambda groups
        file.write('\n')

        addParam('m_lambda', 10.0)                  # mass of l-particles
        addParam('T_lambda', 300)                   # ref. temp. of l-particles
        addParam('tau', 0.1)                        # time constant for thermostat
        addParam('thermostat', 'v-rescale')         # 'v-rescale' or 'langevin'
        addParam('nst_lambda', 100)                 # numSteps between output

        addParam('charge_constraint', 'yes')
        addParam('N_buffers', 1)                    # number of collective buffers

        addParam('m_buf', 10.0)                     # mass of buffer particles
        addParam('multistate_constraint', 'no')     # NOT RELEVANT FOR NOW
        addParam('n_multigroups', 0)                # NOT RELEVANT FOR NOW
        file.write('\n')

        ################ PART 2 - RESIDUE-TYPE SPECIFIC PARAMETERS #############

        def addRes1(name, n_coeffs, dvdl_coeffs, ref_pka):  # formatting function.
            addParam('residue', name)
            addParam('n_coeffs', n_coeffs)

            file.write("{:21s} = ".format('dvdl_coeffs'))
            for coeff in dvdl_coeffs:
                file.write('%.3f ' % (coeff))
            file.write('\n')

            addParam('ref_pka', ref_pka)
            file.write('\n')

        #   resName  numParams  params for ref. potential   refpKa
        addRes1('GLU', 4, [24.685, -577.05, 137.39, -172.69], 4.25)
        addRes1('ASP', 4, [37.822, -566.01, 117.97, -158.79], 3.65)
        # addRes1('BUF', 4, [2010.3, -2023.2, 249.56, -450.63], 4.25)   # old
        addRes1('BUF', 4, [i * countACID for i in [670.1, -674.4, 83.19, -150.21]], 0.00) # new, but might not be necessary in newer commits.

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

        count = 1
        for residue in self.d_residues:      # loop through all residues

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

        # WRITE BUFFER HEAD
        addParam('name', 'BUF')
        addParam('residue_number', 1)
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

################################################################################

    def energy_minimize(self):
        self.__update("energy_minimize", "running gmx grompp and mdrun for energy minimization...")

        os.system("gmx grompp %s >> builder.log 2>&1" % self.g_EM_gromm)
        os.system("gmx mdrun  %s >> builder.log 2>&1" % self.g_EM_md)

    def energy_tcouple(self):
        self.__update("energy_tcouple", "running gmx grompp and mdrun for temperature coupling...")

        os.system("gmx grompp %s >> builder.log 2>&1" % self.g_NVT_gromm)
        os.system("gmx mdrun  %s >> builder.log 2>&1" % self.g_NVT_md)

    def energy_pcouple(self):
        self.__update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")
        os.system("gmx grompp %s >> builder.log 2>&1" % self.g_NPT_gromm)
        os.system("gmx mdrun  %s >> builder.log 2>&1" % self.g_NPT_md)

################################################################################

    def write_run(self, gmxDefaultPath, gmxPhPath):
        self.__update("write_run", "writing run.sh")
        
        with open("run.sh", "w+") as file:
            file.write("#!/bin/bash\n\n")

            if (self.d_constantpH):
                file.write("# source constant-pH gromacs version\n")
                file.write("source %s/bin/GMXRC\n\n" % gmxPhPath)

            file.write("gmx grompp %s\n\n" % self.g_MD_gromm)
            
            if (self.d_constantpH):
                file.write("gmx mdrun %s -nb cpu\n\n" % self.g_MD_md)
            else:
                file.write("gmx mdrun %s\n\n" % self.g_MD_md)

            file.write("# CONTINUE\n")
            file.write("# gmx convert-tpr -s MD.tpr -o MD.tpr -extend <ps>\n")
            file.write("# gmx mdrun %s -cpi state.cpt -append \n\n" % self.g_MD_md)

            if (self.d_constantpH):
                file.write("# source default gromacs version\n")
                file.write("source %s/bin/GMXRC\n" % gmxDefaultPath)

        os.system("chmod +x run.sh")

    def write_reset(self):
        self.__update("write_reset", "writing reset.sh...")

        with open("reset.sh", "w+") as file:
            file.write("#!/bin/bash\n\n")
            
            file.write("rm -rf \\_\\_py* charmm*\n")
            file.write("rm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.cpt *.dat\n")
            file.write("rm -f \\#*\\#\n")
            file.write("rm -f buffer.pdb %s_*.pdb\n" % self.d_pdbName)
            file.write("rm -f run.sh reset.sh jobscript.sh\n")

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
        writeHead("mail-user", "anton.jansen@scilifelab.org")
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

########################## MISCELLANEOUS FUNCTIONS #############################

def backupFile(fname):
    if os.path.isfile(fname):
        num = 1
        while (True):
            if os.path.isfile("#%s.%s#" % (fname, num)):
                num += 1
                continue

            os.system("mv %s '#%s.%s#'" % (fname, fname, num))
            break
