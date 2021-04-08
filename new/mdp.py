import universe
import utils

def gen_mdp(Type, nsteps=25000, nstxout=0, posres=False):
    # HEAD
    if (Type not in ['EM', 'NVT', 'NPT', 'MD']):
        raise Exception("Unknown .mdp Type specified. Types are: IONS, EM, NVT, NPT, MD.")

    utils.update("gen_mdp", "Type={0}, nsteps={1}, nstxout={2}, posres={3}".format(Type, nsteps, nstxout, posres))

    file = open("{0}.mdp".format(Type), 'w')

    def addTitle(title):
        file.write("\n; {0}\n".format(title.upper()))

    def addParam(name, value, comment=''):
        if (comment == ''):
            file.write("{:20s} = {:13s}\n".format(name, str(value)))
        else:
            file.write("{:20s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

    # POSITION RESTRAIN SECTION
    if Type in ['EM', 'MD']:
        if universe.get('d_constantpH') and universe.get('d_restrainpH') and posres:
            addTitle('Position restrain')
            addParam('define', '-DPOSRES -DPOSRES_BUF', 'Position restraints.')
        elif universe.get('d_constantpH') and universe.get('d_restrainpH'):
            addTitle('Position restrain')
            addParam('define', '-DPOSRES_BUF', 'Position restraints.')

    if (Type in ['NVT', 'NPT']): # position restrain temp and press coupling.
        addTitle('Position restrain')
        if universe.get('d_constantpH') and universe.get('d_restrainpH'):
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
    if Type == 'MD' and universe.get('d_restrainpH'):
        addParam('comm-mode', 'Linear', 'Remove center of mass translation.')
        addParam('comm-grps', 'Protein Non-Protein')

    # OUTPUT CONTROL
    addTitle("Output control")
    addParam('nstxout-compressed', nstxout, 'Write frame every %.3f ps.' % (dt * nstxout))

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
    
    if (universe.get('d_modelFF')[0:5].lower() == "charm"): # if we use a CHARMM force field...
        addParam('rcoulomb', 1.2, 'Berk: CHARMM was calibrated for 1.2 nm.')
        addParam('fourierspacing', 0.14, 'Berk: set this to 0.14 for CHARMM.')
    else: # Default for force fields:
        addParam('rcoulomb', 1.0, 'Coulomb cut-off (nm).')

    # VAN DER WAALS
    addTitle("Van der Waals")
    addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')

    if (universe.get('d_modelFF')[0:5].lower() == "charm"): # if we use a CHARMM force field...
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
        addTitle('Generate velocities for startup')
        addParam('gen_vel', 'yes')

    file.close()

    # PUT RELEVANT PARAMETERS IN UNIVERSE
    if (Type == 'MD'):
        universe.add('d_dt', dt)
        universe.add('d_nsteps', nsteps)
        universe.add('d_nstxout', nstxout)
