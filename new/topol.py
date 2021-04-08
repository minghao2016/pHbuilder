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

def generate(d_modelFF, d_modelWater, d_neutralTermini=False):
    universe.add('d_modelFF', d_modelFF)
    universe.add('d_modelWater', d_modelWater)
    universe.add('d_neutralTermini', d_neutralTermini)

    countACID = protein.countRes("ASP") + protein.countRes("GLU")

    # If constant-pH is on,
    if universe.get('d_constantpH'):
        utils.update("generate", 'constant-pH is turned on...')
        
        # and we have at least one protonatable reside,
        if (countACID > 0):
            utils.update("generate", "detected {0} acidic residue(s):".format(countACID))

            count = 1
            for residue in universe.get('d_residues'):
                if (residue.d_resname in ['ASP', 'GLU']):
                    utils.update("generate", "{:3s} {:<4d}".format(residue.d_resname, count))
                count += 1
            utils.update("generate", "(setting protonation state to true for all of these)")

        else:
            utils.update("generate", "no acidic residues detected, turning off constant-pH...")
            universe.add('d_constantpH', False)
            universe.add('d_restrainpH', False)

    else:
        utils.update("generate", 'constant-pH is turned off...')

    utils.update("generate", "using the %s force field with the %s water model..." % (d_modelFF, d_modelWater))

    # Count how many different chains we have
    listChains = []
    # atomCount  = 1
    for residue in universe.get('d_residues'):
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
        xstr += "\n%d" % universe.get('d_constantpH')
    for _ in range(0, len(listChains)):
        xstr += "\n%d\n%d" % (d_neutralTermini, d_neutralTermini)
    xstr += "\nEOF"

    if (d_neutralTermini):
        utils.update("generate", "setting termini to neutral (NH2 and COOH)...")
    else:
        utils.update("generate", "setting termini to charged (NH3+ and COO-) (gmx default)...")
    
    utils.update("generate", "running pdb2gmx to create topol.top...")

    # Generate topology and protonate (make neutral) all GLU and ASP:
    if (d_neutralTermini):
        os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -asp -glu -ignh -ff {2} -water {3} -ter >> builder.log 2>&1 {4}".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_modelFF, d_modelWater, xstr))
    else: # Still a bug here if we have multiple chains, so have to specify the termini multiple times
        os.system("gmx pdb2gmx -f {0} -o {1}_PR2.pdb -asp -glu -ignh -ff {2} -water {3} >> builder.log 2>&1 {4}".format(universe.get('d_nameList')[-1], universe.get('d_pdbName'), d_modelFF, d_modelWater, xstr))

    # Rebuild topology.
    __rebuild_topol()

    # To update d_residues.
    protein.load("{0}_PR2.pdb".format(universe.get('d_pdbName')))

    # To update d_nameList.
    utils.add_to_nameList("{0}_PR2.pdb".format(universe.get('d_pdbName')))

    # To update index.ndx.
    utils.generate_index()

def __rebuild_topol():
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

            # Write the atoms as an extra comment.
            file.write("; {0} {1} {2} {3}\n".format(atomNameList[0], atomNameList[1], atomNameList[2], atomNameList[3]))

            # Atomcount resets for every separate .itp file.
            count = 0
            for residue in universe.get('d_residues'):
                idxList = []
                for atom in residue.d_atoms:
                    # only increase atomcount when we read the relevant chain
                    if residue.d_chain == letter:
                        count += 1

                        if residue.d_resname == resName and atom in atomNameList:
                            idxList.append(count)

                # dihedrals have always four atoms so only activate if we have a length of four.                
                if len(idxList) == 4:
                    utils.update("restrain_dihedrals", "adding restraints for chain {0} {1}-{2}...".format(residue.d_chain, resName, residue.d_resid))
                    file.write("{:<6d} {:<6d} {:<6d} {:<6d}  {}  {}  {}  {}\n".format(idxList[0], idxList[1], idxList[2], idxList[3], Type, phi, dphi, fc))
