import os, universe, utils, protein

def add_mol(itpfname, comment, molname=None, molcount=None):
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

def generate(d_modelFF, d_modelWater, d_terministring="", d_customitp=""):
    # Internal helper function.
    def rebuild_topol():
        # If we have only one chain, gromacs will put everything in topol.top.
        # If we have more than one chain, gromacs will do it for us.
        if (len(universe.get('d_chain')) <= 1):
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
            file.write("#include \"{0}.ff/forcefield.itp\"\n\n".format(universe.get('d_modelFF')))

            file.write("; Include protein topology\n")
            for letter in universe.get('d_chain'):
                file.write("#include \"topol_Protein_chain_{0}.itp\"\n".format(letter))
            file.write('\n')

            file.write('[ system ]\n')
            file.write('{0}\n\n'.format(universe.get('d_pdbName')))

            file.write('[ molecules ]\n')
            file.write('; Compounts \t\t #mols\n')
            for letter in universe.get('d_chain'):
                file.write("Protein_chain_{0}\t\t1\n".format(letter))

    # ADD RELEVANT PARAMETERS TO UNIVERSE ######################################
    universe.add('d_modelFF', d_modelFF)
    universe.add('d_modelWater', d_modelWater)
    universe.add('d_terministring', d_terministring)

    # USER UPDATE STUFF ########################################################
    countACID = protein.countRes("ASP") + protein.countRes("GLU")

    # If constant-pH is on,
    if universe.get('ph_constantpH'):
        utils.update("generate", 'constant-pH is turned on...')
        
        # and we have at least one protonatable reside,
        if countACID > 0:
            utils.update("generate", "detected {} acidic residue(s):".format(countACID))

            for letter in universe.get('d_chain'):
                count = 1

                for residue in universe.get('d_residues'):
                    if residue.d_chain == letter:
                        count += 1

                        if residue.d_resname in ['ASP', 'GLU']:
                            utils.update("generate", "{:3s}-{:<4d} in chain {}".format(residue.d_resname, count, letter))

            utils.update("generate", "(setting protonation state to true (option 1) for all of these)")

        else:
            utils.update("generate", "no acidic residues detected, constant-pH is turned off...")
            universe.add('ph_constantpH', False)
            universe.add('ph_restrainpH', False)

    else:
        utils.update("generate", 'constant-pH is turned off...')
        universe.add('ph_restrainpH', False) # If ph_constantpH is False then this is also False.

    utils.update("generate", "using the {} force field with the {} water model...".format(d_modelFF, d_modelWater))

    # CREATE EOFSTRING FOR PDB2GMX COMMAND #####################################

    # Here we create EOFstring to circumvent pd2gmx prompting for user input.
    xstr = "<< EOF"
    # The for-loop sets the protonation state of all GLUs and ASPs to 1
    # (True) if ph_constantpH is True, and to 0 (False) if ph_constantpH is False.
    for chain in universe.get('d_chain'):
        for residue in universe.get('d_residues'):
            if residue.d_resname in ['ASP', 'GLU'] and residue.d_chain == chain:
                xstr += "\n{:d}".format(universe.get('ph_constantpH'))
        
        # Furthermore, if the user specified a two-letter string for the termini,
        # add those to the EOFstring as well:
        if (d_terministring != ""):
            xstr += "\n{}".format(d_terministring[0])
            xstr += "\n{}".format(d_terministring[1])
    # End EOFstring:
    xstr += "\nEOF"
    # print(xstr) # Debug

    # UPDATE USER ABOUT WHAT WE DO WITH THE TERMINI ############################

    if (d_terministring != ""):
        utils.update("generate", "Using options {} for termini...".format(d_terministring))
    else:
        utils.update("generate", "No termini specified, using gmx default (00 = NH3+ and COO-)...")
    
    utils.update("generate", "running pdb2gmx to create {}_PR2.pdb and topol.top...".format(universe.get('d_pdbName')))

    # RUN ACTUAL PDB2GMX COMMAND ###############################################

    if (d_terministring != ""):
        os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -asp -glu -ignh -ff {2} -water {3} -ter >> builder.log 2>&1 {4}".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_modelFF, d_modelWater, xstr))
    else:
        os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -asp -glu -ignh -ff {2} -water {3} >> builder.log 2>&1 {4}".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_modelFF, d_modelWater, xstr))

    # Rebuild topology.
    rebuild_topol()

    # Use custom topol_Protein_chain_A.itp (this is a temporary fix for charmm36-mar2019-m4)
    if (d_customitp != ""):
        utils.update("generate", "Overwriting topol_Protein_chain_A.itp generated by pdb2gmx with {}...".format(d_customitp))
        os.system("cp {} .".format(d_customitp))

    # To update d_residues.
    protein.load("{0}_PR2.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_PR2.pdb".format(universe.get('d_pdbName')))

def restrain_dihedrals(resName, atomNameList, Type, phi, dphi, fc):
    utils.update("restrain_dihedrals", "will add restraints for {0} (all chains)...".format(resName))

    # Every chain has its own .itp file, so we loop through every file:
    for letter in universe.get('d_chain'):
        # This is to make sure we don't have multiple headers when we add multiple different restraints.
        first = False
        if not "[ dihedral_restraints ]" in open("topol_Protein_chain_{0}.itp".format(letter)).read():
            first = True

        # Append to the end of relevant .itp file:
        with open("topol_Protein_chain_{0}.itp".format(letter), 'a') as file:
            if first:
                file.write("[ dihedral_restraints ]\n")
                file.write("; ai aj ak al type phi dphi fc\n")

            # Atomcount resets for every separate .itp file.
            atomCount = 0
            for residue in universe.get('d_residues'):
                # Dictionary resets for every residue, as we can of course have
                # multiple ASPs or GLUs in one chain.
                dictionary = {}
                for atom in residue.d_atoms:
                    # Only increase atomcount when we read the relevant chain:
                    if residue.d_chain == letter:
                        atomCount += 1

                        if residue.d_resname == resName and atom in atomNameList:
                            dictionary[atom] = atomCount

                if len(dictionary) == 4:
                    utils.update("restrain_dihedrals", "adding restraints for chain {0} {1}-{2}...".format(residue.d_chain, resName, residue.d_resid))
                    for atom in atomNameList:
                        file.write("{:<6d} ".format(dictionary[atom]))
                    file.write(" {}  {}  {}  {}\n".format(Type, phi, dphi, fc))

def restrain_dihedrals_by_idx(indices, Type, phi, dphi, fc):
    utils.update("restrain_dihedrals", "will restrain dihedral {} {} {} {} (Type={}, phi={}, dphi={}, fc={})".format(indices[0], indices[1], indices[2], indices[3], Type, phi, dphi, fc))

    # Every chain has its own .itp file, so we loop through every file:
    for letter in universe.get('d_chain'):
        # This is to make sure we don't have multiple headers when we add multiple different restraints.
        first = False
        if not "[ dihedral_restraints ]" in open("topol_Protein_chain_{0}.itp".format(letter)).read():
            first = True

        # Append to the end of relevant .itp file:
        with open("topol_Protein_chain_{0}.itp".format(letter), 'a') as file:
            if first:
                file.write("[ dihedral_restraints ]\n")
                file.write("; ai aj ak al type phi dphi fc\n")

            for idx in indices:
                file.write("{:<6d} ".format(idx))
                
            file.write(" {}  {}  {}  {}\n".format(Type, phi, dphi, fc))
