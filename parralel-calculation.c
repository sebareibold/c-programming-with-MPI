#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>
#include <mpi.h>

#define Cx 0.1f
#define Cy 0.1f

double sampleTime()
{
    struct timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return ((double)tv.tv_sec + ((double)tv.tv_nsec) / 1000000000.0);
}

int main(int argc, char *argv[])
{

    MPI__Init();

    MPI_Finalize();
}