#!bin/bash
mpicc -g blumecapel_pt_mpi.c -lm -O3 -std=c99
for((c=0;c<10;c=c+1))
do
    mpirun -np $1 ./a.out $2 $3 $4 $5 $6 $7 $c
done;