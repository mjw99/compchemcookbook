export AMBERHOME=/netscratch/mw529/code/AMBER/amber11

#export MPI_RUN:=/server-home/netbin/mpi/mpich2-1.0.7-ifort-10.1.008/bin/mpirun -np 2

#The actual MD commands
export SANDER:=$(AMBERHOME)/bin/sander -O
export PMEMD:=$(AMBERHOME)/bin/pmemd.cuda -O
#export S:=$(MPI_RUN) $(AMBERHOME)/exe/pmemd.MPI -O  < /dev/null 

#Macro to genearate a pdb
export MAKEPDB:=$(AMBERHOME)/bin/ambpdb -p prmtop < restrt > snap.pdb

