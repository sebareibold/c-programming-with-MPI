#!/bin/bash

# Definimos el número de procesos 
NUM_PROCESOS=7

#Compilamos el programa MPI
mpiCC parralel-calculation-v-s-c.c -o parralel-calculation

# Ejecuta el programa MPI 
mpirun -np $NUM_PROCESOS ./parralel-calculation 100  3
