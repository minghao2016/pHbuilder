import os, universe, utils

from mdp import gen_mdp
from constantph import gen_constantpH

def energy_minimize():
    gen_mdp('EM')
    
    utils.update("energy_minimize", "running gmx grompp and mdrun for energy minimization...")

    os.system("gmx grompp -f EM.mdp -c {0} -p topol.top -o EM.tpr -r {0} >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx mdrun -deffnm EM -c {0}_EM.pdb >> builder.log 2>&1".format(universe.get('d_pdbName')))

    utils.add_to_nameList("{0}_EM.pdb".format(universe.get('d_pdbName')))

def energy_tcouple():
    gen_mdp('NVT')
    
    utils.update("energy_tcouple", "running gmx grompp and mdrun for temperature coupling...")

    os.system("gmx grompp -f NVT.mdp -c {0} -p topol.top -o NVT.tpr -r {0} >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx mdrun deffnm NVT -c {0}_NVT.pdb >> builder.log 2>&1".format(universe.get('d_pdbName')))

    utils.add_to_nameList("{0}_NVT.pdb".format(universe.get('d_pdbName')))

def energy_pcouple():
    gen_mdp('NPT')

    utils.update("energy_pcouple", "running gmx grompp and mdrun for pressure coupling...")

    os.system("gmx grompp -f NPT.mdp -c {0} -p topol.top -o NPT.tpr -r {0} >> builder.log 2>&1".format(universe.get('d_nameList')[-1]))
    os.system("gmx mdrun -deffnm NPT -c {0}_NPT.pdb >> builder.log 2>&1".format(universe.get('d_pdbName')))    

    utils.add_to_nameList("{0}_NPT.pdb".format(universe.get('d_pdbName')))
