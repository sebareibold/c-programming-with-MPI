#!/bin/bash

# Definimos el número de procesos 
NUM_PROCESOS=4

#Compilamos el programa MPI
mpiCC parralel-calculation.c -o parralel-calculation


# Ejecuta el programa MPI 
mpirun -np $NUM_PROCESOS ./parralel-calculation 16  1