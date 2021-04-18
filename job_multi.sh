#!/bin/bash
#SBATCH --time=1-23:59:59
#SBATCH --nodes=1
#SBATCH --partition=lindahl
#SBATCH --job-name=ASP_tri
#SBATCH --mail-user=anton.jansen@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -C gpu --gres=gpu:4

# Load modules:
module load cmake/latest
module load gcc/7.4
module load cuda/10.2

# Compile our custom Gromacs version on cluster backend node:
mkdir build
cd build
CC=gcc-7 CXX=g++-7 cmake ~/gromacs-constantph -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=${PWD}/.. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA
make -j 12
make install -j 12
cd ..
rm -r build
source ${PWD}/bin/GMXRC

# Generate MD.tpr:
gmx grompp -f MD.mdp -c ASP_tri_NPT.pdb -p topol.top -n index.ndx -o MD.tpr -r ASP_tri_NPT.pdb

# Start four runs:
cpus=$(( $SLURM_JOB_CPUS_PER_NODE / 4 ))
for i in 0 1 2 3; do
    mkdir run${i}
    cp MD.tpr run${i}
    cd run${i}
    export CUDA_VISIBLE_DEVICES=$i
    gmx mdrun -deffnm MD -x MD.xtc -c ASP_tri_MD.pdb -pme cpu -nt $cpus -pin on -pinstride 1 -pinoffset $(( $cpus * $i )) &
    cd ..
done
wait
