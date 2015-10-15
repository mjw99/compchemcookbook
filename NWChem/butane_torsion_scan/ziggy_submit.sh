#!/bin/bash
# Example submit script to use Ziggy nwchem


NAME=IJDNQMDRQITEOD-UHFFFAOYSA-N

/opt/torque/bin/qsub << EOF

########################
# Start PBS directives #
########################

#          Set the name of the job (up to 15 characters, 
#          no blank spaces, start with alphanumeric character)

#PBS -N $NAME


#          Specify the queue. 

#PBS -q s16


#          Specify the maximum wall clock time.
#          Specify the number of nodes requested and the
#          number of processors per node. 


#PBS -l nodes=1:ppn=16,walltime=00:30:00

#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output. 

#PBS -j oe

########################
# End PBS directives   #
########################





module purge
module add gcc/4.8.3
module add mpi/openmpi/gnu/1.8.1
module add nwchem


# Change to directory where you submitted the job from
cd \${PBS_O_WORKDIR}

# Needed so that python picks up the scripts in this directory
export PYTHONPATH=\${PYTHONPATH}:\${PBS_O_WORKDIR}

# Run the job
mpiexec -np 16 nwchem $NAME.nwin > $NAME.nwlog

# Clean ups
rm -f nwchem.py*
rm -f *.hess
rm -f *.movecs
rm -f *.db

EOF
