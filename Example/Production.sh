#!/bin/bash

# Change this path to the one you had installed Gromacs.
source /usr/local/gromacs/bin/GMXRC

gmx_mpi grompp -f production.mdp -o Example.tpr -c start.gro -p system.top -n index.ndx -maxwarn 1
gmx_mpi mdrun -v -deffnm Example -plumed plumed.dat
