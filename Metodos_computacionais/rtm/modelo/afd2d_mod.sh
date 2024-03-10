#!/bin/bash
parfile=parfile_rtm.txt
mpirun -n 4 afd_mpi.x<<EOF
${parfile}
EOF
