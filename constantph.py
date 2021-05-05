import os, universe, utils, md

def gen_constantpH(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE, cal=False, lambdaInit=0.5):
    # Hardcoded stuff
    GLU_pKa   = 4.25
    GLU_atoms = [' CG ', ' CD ', ' OE1', ' OE2', ' HE2'] # atoms part of model
    GLU_qqA   = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    GLU_qqB   = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge
    
    ASP_pKa   = 3.65
    ASP_atoms = [' CB ', ' CG ', ' OD1', ' OD2', ' HD2'] # atoms part of model
    ASP_qqA   = [-0.21 ,  0.75 ,  -0.55,  -0.61,  0.44 ] # protonated charge
    ASP_qqB   = [-0.28 ,  0.62 ,  -0.76,  -0.76,  0.00 ] # deprotonated charge

    BUF_qqA   = [-0.0656, 0.5328, 0.5328]
    BUF_qqB   = [-0.8476, 0.4238, 0.4238]

    # Skip this entire step if ph_constantpH is false.
    if (not universe.get('ph_constantpH')):
        utils.update("generate_phdata", "skipping this step...")
        return

    # Load dV/dl coefficients
    GLU_dvdl = universe.get('ph_GLU_dvdl')
    ASP_dvdl = universe.get('ph_ASP_dvdl')

    if (universe.get('ph_restrainpH')):
        BUF_dvdl = universe.get('ph_BUF_dvdl')

    # Check if MD.mdp exists.
    if (not os.path.isfile("MD.mdp")):
        utils.update("generate_phdata", "MD.mdp does not exist, creating...")
        md.gen_mdp('MD', universe.get('d_nsteps'), universe.get('d_nstxout'))

    # Generate default index file.
    utils.generate_index()

    # Update user.
    utils.update("generate_phdata", "ph_pH={}, ph_lambdaM={}, ph_nstout={}, ph_barrierE={}".format(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE))

    file = open('MD.mdp', 'a')

    # Formatting function.
    def addParam(name, value, comment = "NUL"):
        if (comment == "NUL"):
            file.write("{:54s} = {:13s}\n".format(name, str(value)))
        else:            
            file.write("{:54s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    file.write("\n; CONSTANT PH\n")

    # PART 1 - WRITE GENERAL PARAMETERS ########################################
    
    addParam('lambda-dynamics', 'yes')
    addParam('lambda-dynamics-simulation-ph', ph_pH)
    addParam('lambda-dynamics-lambda-particle-mass', ph_lambdaM)
    addParam('lambda-dynamics-update-nst', ph_nstout)
    addParam('lambda-dynamics-tau', 2.0) # hardcoded

    if cal:
        addParam('lambda-dynamics-calibration', 'yes')

    if universe.get('ph_restrainpH'):
        addParam('lambda-dynamics-charge-constraints', 'yes')

    # Compile a list of acidic residues and their ResIDs.
    acidicResidueNameList = []; acidicResidueNumberList = []
    acidicResidueTypeList = []
    
    for residue in universe.get('d_residues'):
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

    if universe.get('ph_restrainpH'):         # If we restrain the charge 
        acidicResidueTypeList.append('BUF')   # we also have BUF.

    addParam('lambda-dynamics-number-lambda-residues', len(acidicResidueTypeList))
    
    if universe.get('ph_restrainpH'):
        addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList) + 1)
    else:
        addParam('lambda-dynamics-number-atom-collections', len(acidicResidueNameList))

    file.write('\n')

    # print(acidicResidueNameList)   # debug
    # print(acidicResidueNumberList) # debug
    # print(acidicResidueTypeList)   # debug

    # PART 2 - WRITE RESIDUE-TYPE SPECIFIC STUFF ###############################

    def writeBlock(number, name, dvdl, pKa, ph_barrierE, qqA, qqB):

        def to_string(Input):
            string = ""
            for element in Input:
                string += str(element)
                string += " "
            return string

        addParam('lambda-dynamics-residue%s-name'              % (number), name)
        addParam('lambda-dynamics-residue%s-dvdl-coefficients' % (number), to_string(dvdl))
        addParam('lambda-dynamics-residue%s-reference-pka'     % (number), pKa)
        addParam('lambda-dynamics-residue%s-barrier'           % (number), ph_barrierE)
        addParam('lambda-dynamics-residue%s-charges-state-A'   % (number), to_string(qqA))
        addParam('lambda-dynamics-residue%s-charges-state-B'   % (number), to_string(qqB))
        
        file.write('\n')

    for idx in range(0, len(acidicResidueTypeList)):
        if (acidicResidueTypeList[idx] == 'GLU'):
            writeBlock(idx + 1, 'GLU', GLU_dvdl, GLU_pKa, ph_barrierE, GLU_qqA, GLU_qqB)

        if (acidicResidueTypeList[idx] == 'ASP'):
            writeBlock(idx + 1, 'ASP', ASP_dvdl, ASP_pKa, ph_barrierE, ASP_qqA, ASP_qqB)

        if (acidicResidueTypeList[idx] == 'BUF'):
            # Multiplication is no-longer necessary because of Paul's commit on January 25th:
            # writeBlock(idx + 1, 'BUF', [i * len(acidicResidueNameList) for i in BUF_dvdl], 0, 0, BUF_qqA, BUF_qqB)
            writeBlock(idx + 1, 'BUF', BUF_dvdl, 0, 0, BUF_qqA, BUF_qqB)

    # PART 3 - WRITE INDIVIDUAL RESIDUE/LAMBDA-GROUP STUF ######################

    def writeResBlock(number, name, indexLambda, indexName):
        addParam('lambda-dynamics-atom-set%s-name'                  % (number), name)
        addParam('lambda-dynamics-atom-set%s-lambda-residues-index' % (number), indexLambda)
        addParam('lambda-dynamics-atom-set%s-index-group-name'      % (number), indexName)
        addParam('lambda-dynamics-atom-set%s-initial-lambda'        % (number), lambdaInit)
        
        if universe.get('ph_restrainpH'):
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

    if universe.get('ph_restrainpH'):
        writeResBlock(
                        len(acidicResidueNameList) + 1,
                        'BUF',
                        acidicResidueTypeList.index('BUF') + 1,
                        'LAMBDA%s' % (len(acidicResidueNameList) + 1)
                        )

    file.close() # MD.mdp

    # PART 4 - APPEND THE LAMBDA INDEX GROUPS TO INDEX.NDX #####################

    file = open('index.ndx', 'a') # Append to existing index.ndx

    # Function for adding an indexList to index.ndx
    def writeTheGroup(number, indexList):
        file.write('\n[ LAMBDA{} ]\n'.format(number))
        for index in indexList:
            file.write('{} '.format(index))
        file.write('\n')

    atomCount = 1   # Keeps track of the atom number.
    grpNum    = 1   # Keeps track of the group (the LAMBDA%s).
    
    for residue in universe.get('d_residues'):       # loop through all residues

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

    if universe.get('ph_restrainpH'):

        atomCount = 1; indexList = []

        for residue in universe.get('d_residues'):
            for atom in residue.d_atoms:
            
                if (residue.d_resname == 'BUF'):
                    indexList.append(atomCount)
            
                atomCount += 1
        
        writeTheGroup(grpNum, indexList)

    file.close() # index.ndx

    # Put relevant pH variables in universe
    universe.add('ph_pH', ph_pH)
    universe.add('ph_lambdaM', ph_lambdaM)
    universe.add('ph_nstout', ph_nstout)
    universe.add('ph_barrierE', ph_barrierE)
