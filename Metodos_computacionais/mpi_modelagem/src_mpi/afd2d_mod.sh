#!/bin/bash
parfile=afd2d_parfile.txt
mpirun -n 4 afd_mpi.x<<EOF
${parfile}
EOF
