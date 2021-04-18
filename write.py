import os
import utils
import universe

def reset():
    utils.update("reset", "writing reset.sh...")

    with open("reset.sh", "w+") as file:
        file.write("#!/bin/bash\n\n")
        
        file.write("if [ -f \"%s_MD.pdb\" ]\nthen\n" % universe.get('d_pdbName'))
        file.write("\tread -p \"Warning: simulation has finished. Proceed? (y)\" var\n")
        file.write("else\n")
        file.write("\trm -rf \\_\\_py* charmm*\n")
        file.write("\trm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg\n")
        file.write("\trm -f \\#*\\#\n")
        file.write("\trm -f buffer.pdb %s_*.pdb\n" % universe.get('d_pdbName'))
        file.write("\trm -f run.sh reset.sh jobscript.sh universe\n")                   
        file.write("fi\n\n")

        file.write("if [ \"${var}\" = \"y\" ]\nthen\n")
        file.write("\trm -rf \\_\\_py* charmm*\n")
        file.write("\trm -f *.itp *.top *.mdp *.tpr *.log *.ndx *.edr *.trr *.xtc *.cpt *.dat *.pdf *.xvg\n")
        file.write("\trm -f \\#*\\#\n")
        file.write("\trm -f buffer.pdb %s_*.pdb\n" % universe.get('d_pdbName'))
        file.write("\trm -f run.sh reset.sh jobscript.sh universe\n")            
        file.write("fi\n\n")

    os.system("chmod +x reset.sh")

def run(gmxPath="/usr/local/gromacs", options=""):
    utils.update("run", "gmxPath={0}, additional options= {1}".format(gmxPath, options))
    
    with open("run.sh", 'w') as file:
        file.write("#!/bin/bash\n\n")

        file.write("# Gromacs version to use:\n")
        file.write("source {0}/bin/GMXRC\n\n".format(gmxPath))

        file.write("gmx grompp -f MD.mdp -c {0} -p topol.top -n index.ndx -o MD.tpr -r {0}\n".format(universe.get('d_nameList')[-1]))
        file.write("gmx mdrun -deffnm MD -v -c {0}_MD.pdb -x MD.xtc {1}\n".format(universe.get('d_pdbName'), options))

    os.system("chmod +x run.sh")

def jobscript(jobName, jobTime, nodes, ntasks, queue):
    utils.update("jobscript", "jobName={0}, jobTime={1}(hrs), nodes={2}, ntasks={3}, queue={4}...".format(jobName, jobTime, nodes, ntasks, queue))

    file = open("jobscript.sh", 'w')
    
    def writeHead(param, value):
        file.write("#SBATCH --%s=%s\n" % (param, value))

    def moduleLoad(value):
        file.write("module load {0}\n".format(value))

    file.write("#!/bin/bash\n")

    writeHead("time", "%d-%.2d:00:00" % (int(jobTime / 24), jobTime % 24))
    writeHead("nodes", nodes)
    writeHead("ntasks", ntasks)
    writeHead("partition", queue)
    writeHead("job-name", jobName)
    writeHead("mail-user", "anton.jansen@scilifelab.se")
    writeHead("mail-type", "ALL")
    file.write("#SBATCH -C gpu --gres=gpu:4\n")

    moduleLoad("cmake/latest")
    moduleLoad("gcc/7.4")
    moduleLoad("cuda/10.2")

    file.write('\n')

    if universe.get('ph_constantpH'):
        file.write("\n# compile our custom Gromacs version on cluster backend node\n")
        file.write("mkdir build\n")
        file.write("cd build\n")
        file.write("CC=gcc-7 CXX=g++-7 cmake ~/gromacs-constantph -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=${PWD}/.. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA\n")
        file.write("make -j 12\n")
        file.write("make install -j 12\n")
        file.write("cd ..\n")
        file.write("rm -r build\n")
        file.write("source ${PWD}/bin/GMXRC\n\n")
    else:
        file.write("module load gromacs/2021.1\n\n")

    file.write("gmx grompp -f MD.mdp -c {0} -p topol.top -n index.ndx -o MD.tpr -r {0}\n".format(universe.get('d_nameList')[-1]))
    
    if universe.get('ph_constantpH'):
        file.write("gmx mdrun -deffnm MD -c {0}_MD.pdb -x MD.xtc -pme cpu\n".format(universe.get('d_pdbName')))
    else:
        file.write("gmx mdrun -deffnm MD -c {0}_MD.pdb -x MD.xtc\n".format(universe.get('d_pdbName')))
