#!/bin/bash
#
# LSF batch script to run the test MPI code
#
#PBS -N HRLDAS-restart-test
#PBS -A UWIS0037                                     
#PBS -l walltime=04:00:00                            
#PBS -q regular                                      
#PBS -j oe                                             
#PBS -l select=1:ncpus=36:mpiprocs=12
#PBS -m abe                                            
#PBS -M gaalexander3@wisc.edu			
#
cd /glade/u/home/galexand/work/noahmp_hue/hrldas/run
mpiexec_mpt ./hrldas.exe >& milwaukee-gfortran-mpi.log


