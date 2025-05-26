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
    register int p, i, j;
    if (argc != 3)
    {
        printf("Parámetros: Tlado pasos\n");
        exit(1);
    }

    int Tlado = atoi(argv[1]); // N --> Tamaño de la matriz cuadrada
    int pasos = atoi(argv[2]); // M --> Cantidad de pasos

    if (Tlado <= 0 || pasos < 0)
    {
        printf("ERROR: valores incorrectos en parámetros\n");
        exit(1);
    }

    int columnas = Tlado;
    int filas = Tlado;

    // Reservamos la memoria de las filas y columnas
    int *fila_arriba = malloc(sizeof(float) * columnas);
    int *fila_abajo = malloc(sizeof(float) * columnas);
    int *columna_izquierda = malloc(sizeof(float) * filas);
    int *columna_derecha = malloc(sizeof(float) * filas);

    // Seteamos dichas columnas y filas con valores 0
    bzero(fila_arriba, sizeof(float) * columnas);
    bzero(fila_abajo, sizeof(float) * columnas);
    bzero(columna_izquierda, sizeof(float) * filas);
    bzero(columna_derecha, sizeof(float) * filas);

    float **matrizLocal;
    float **matrizSiguiente;
    float **aux; // Puntero auxiliar para intercambio de punteros

    // Reserva de espacio para la matriz Actual
    matrizLocal = (float **)malloc(sizeof(float *) * Tlado);
    *matrizLocal = (float *)malloc(sizeof(float) * Tlado * Tlado);

    // Reserva de espacio para la matriz Siguiente
    matrizSiguiente = (float **)malloc(sizeof(float *) * Tlado);
    *matrizSiguiente = (float *)malloc(sizeof(float) * Tlado * Tlado);

    if (matrizLocal == NULL || matrizSiguiente == NULL)
    {
        printf("ERROR: No se pudo reservar memoria\n");
        exit(2);
    }

    // Asigna las direcciones de inicio para cada fila de las matrices matrizLocal y matrizSiguiente
    // (desde la fila 1 en adelante) dentro del único bloque contiguo de memoria asignado previamente.
    // Forma eficiente de manejar matrices 2D en C. (Arma la matris)
    for (i = 1; i < Tlado; i++)
    {
        matrizLocal[i] = *matrizLocal + Tlado * i;
        matrizSiguiente[i] = *matrizSiguiente + Tlado * i;
    }

    // Inicialización de la matriz Actual
    for (i = 0; i < Tlado; i++)
        for (j = 0; j < Tlado; j++)
            matrizLocal[i][j] = (float)((i + 1)) * (Tlado + i) * (j + 1) * (Tlado + j);

    // Ejecución de los pasos de simulación
    double time_spent = sampleTime();

    // ARRANCAMOS a pensar en MPI

    int rank, size, rank_cart,
        arriba, abajo, izq, der; // Vecinos en la topologia 2D.

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm COMM_CART;
    MPI_Status status;

    // Damos la cantidad de dimensiones (columnas y filas)
    int dims[] = {((int)sqrt(size)), size / ((int)sqrt(size))};
    int periods[] = {0, 0}; // No sera periodica ninguna dimension.
    float subMatriz();
    int coord[2]; // Podemos almacenar coordenadas 2D

    int tamañio_local_fil = Tlado / dims[0];
    int tamañio_local_col = Tlado / dims[1];

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &COMM_CART);
    if (COMM_CART != MPI_COMM_NULL)
    {
        MPI_Comm_rank(COMM_CART, &rank_cart); // Obtenemos el nuevo rango, por el nuevo comunicador

        MPI_Type_vector(tamañio_local_fil, 1, tamañio_local_col, MPI_INT, &subMatriz); // Guardamos el tipo de dato
        MPI_Type_commit(&subMatriz);                                                   // Confirmamos el tipo para poder usarlo.

        if (rank == 0)
        {

            int aux, rank_cart;
            for (int i = 0; i < dims[0]; i++) // para cada proceso se reparte su parte
                for (int j = 0; j < dims[1]; j++)
                {
                    coord[0] = i;
                    coord[1] = j;
                    MPI_Cart_rank(COMM_CART, coord, &aux);

                    // Sirve para mandar a todos menos al procs 0
                    if (aux != rank_cart)
                    {
                        // enviamos las submatriz
                        MPI_Send(&matrizLocal[i * tamañio_local_fil][j * tamañio_local_col], 1, subMatriz, aux, 0, COMM_CART);
                    }
                }

            MPI_Cart_coords(COMM_CART, rank_cart, 2, coord);
        }
        else if (rank != 0)
        {
            MPI_Cart_coords(COMM_CART, rank_cart, 2, coord);
            MPI_Recv(&matrizLocal[tamañio_local_fil * coord[0]][tamañio_local_col * coord[1]], 1, subMatriz, MPI_ANY_SOURCE, 0, COMM_CART, &status);
        }
    }

    // DESPUES DE CADA PODEMOS PENSARLO EN MANDAR LAS COLUMNAS SI ES Q TIENE VECINOS, COMO EL EJEMPLO DEL UNIVERSIDAD DE GRANADA
    MPI_Datatype vectorVertical;
    MPI_Type_vector(Tlado, 1, 1, MPI_INT, &vectorVertical);
    MPI_Type_commit(&vectorVertical);

    MPI_Cart_shift(COMM_CART, 0, 1, &arriba, &abajo); // dame los vecinos (destino)
    // lo mismo para la dimension 1
    MPI_Cart_shift(COMM_CART, 1, 1, &izq, &der);

    if (arriba != MPI_PROC_NULL) // si tengo vecino arriba mando mi fila superior
        MPI_Send(&ImagenLocal[tamañio_local_fil * coord[0]]
                             [tamañio_local_col * coord[1]],
                 tamañio_local_fil, MPI_INT, arriba, 0, COMM_CART);

    if (abajo != MPI_PROC_NULL) // si tengo vecino abajo recibo su fila superior
        MPI_Recv(&ImagenLocal[(tamañio_local_fil * coord[0]) + tamañio_local_fil]
                             [tamañio_local_col * coord[1]],
                 tamañio_local_fil, MPI_INT, abajo, 0, COMM_CART, &status);

    if (abajo != MPI_PROC_NULL) // si tengo vecino abajo mando mi fila inferior
        MPI_Send(&ImagenLocal[(tamañio_local_fil * coord[0]) + tamañio_local_fil - 1]
                             [tamañio_local_col * coord[1]],
                 tamañio_local_fil, MPI_INT, abajo, 0, COMM_CART);

    if (arriba != MPI_PROC_NULL) // si tengo vecino arriba recibo su fila inferior
        MPI_Recv(&ImagenLocal[(tamañio_local_fil * coord[0]) - 1]
                             [tamañio_local_col * coord[1]],
                 tamañio_local_fil, MPI_INT, arriba, 0, COMM_CART, &status);

    if (izq != MPI_PROC_NULL) // si tengo vecino izquierda mando mi columna izquida.
        MPI_Send(&ImagenLocal[tamañio_local_fil * coord[0]]
                             [tamañio_local_col * coord[1]],
                 1, vectorVertical, izq, 0, COMM_CART);

    if (der != MPI_PROC_NULL) // si tengo vecino derecha recibo su columna derecha.
        MPI_Recv(&ImagenLocal[tamañio_local_fil * coord[0]]
                             [(tamañio_local_col * coord[1]) + tamañio_local_col],
                 1, vectorVertical, der, 0, COMM_CART,
                 &status);

    if (der != MPI_PROC_NULL) // si tengo vecino derecha mando mi columna derecha
        MPI_Send(&ImagenLocal[tamañio_local_fil * coord[0]]
                             [(tamañio_local_col * coord[1]) + tama_local - 1],
                 1, vectorVertical, der, 0,
                 COMM_CART);

    if (izq != MPI_PROC_NULL) // si tengo vecino izquierda recibo su columna derecha.
        MPI_Recv(&ImagenLocal[tamañio_local_fil * coord[0]]
                             [(tamañio_local_col * coord[1]) - 1],
                 1, vectorVertical, izq, 0, COMM_CART, &status);

    // REALIZAR LOS CALCULOS

    float e_arriba, e_abajo, e_izq, e_der, yo;
    for (p = 0; p < pasos; p++)
    {

        // Se procesa el interior
        for (i = 1; i < filas - 1; i++)
            for (j = 1; j < columnas - 1; j++)
            {
                yo = matrizLocal[i][j];
                e_arriba = matrizLocal[i - 1][j];
                e_abajo = matrizLocal[i + 1][j];
                e_izq = matrizLocal[i][j - 1];
                e_der = matrizLocal[i][j + 1];
                matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
            }

        // Se procesa fila superior
        i = 0;
        for (j = 1; j < columnas - 1; j++)
        {
            yo = matrizLocal[i][j];
            e_arriba = fila_arriba[j];
            e_abajo = matrizLocal[i + 1][j];
            e_izq = matrizLocal[i][j - 1];
            e_der = matrizLocal[i][j + 1];
            matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
        }

        // Se procesa fila inferior
        i = filas - 1;
        for (j = 1; j < columnas - 1; j++)
        {
            yo = matrizLocal[i][j];
            e_arriba = matrizLocal[i - 1][j];
            e_abajo = fila_abajo[j];
            e_izq = matrizLocal[i][j - 1];
            e_der = matrizLocal[i][j + 1];
            matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
        }

        // Se procesa columna izquierda
        j = 0;
        for (i = 1; i < filas - 1; i++)
        {
            yo = matrizLocal[i][j];
            e_arriba = matrizLocal[i - 1][j];
            e_abajo = matrizLocal[i + 1][j];
            e_izq = columna_izquierda[i];
            e_der = matrizLocal[i][j + 1];
            matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
        }

        // Se procesa columna derecha
        j = columnas - 1;
        for (i = 1; i < filas - 1; i++)
        {
            yo = matrizLocal[i][j];
            e_arriba = matrizLocal[i - 1][j];
            e_abajo = matrizLocal[i + 1][j];
            e_izq = matrizLocal[i][j - 1];
            ;
            e_der = columna_derecha[i];
            matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
        }

        // Se procesa esquina superior izquierda
        i = 0;
        j = 0;
        yo = matrizLocal[i][j];
        e_arriba = fila_arriba[0];
        e_abajo = matrizLocal[i + 1][j];
        e_izq = columna_izquierda[0];
        e_der = matrizLocal[i][j + 1];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);

        // Se procesa esquina superior derecha
        i = 0;
        j = columnas - 1;
        yo = matrizLocal[i][j];
        e_arriba = fila_arriba[columnas - 1];
        e_abajo = matrizLocal[i + 1][j];
        e_izq = matrizLocal[i][j - 1];
        e_der = columna_derecha[i];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);

        // Se procesa esquina inferior izquierda
        i = filas - 1;
        j = 0;
        yo = matrizLocal[i][j];
        e_arriba = matrizLocal[i - 1][j];
        e_abajo = fila_abajo[j];
        e_izq = columna_izquierda[filas - 1];
        e_der = matrizLocal[i][j + 1];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);

        // Se procesa esquina inferior derecha
        i = filas - 1;
        j = columnas - 1;
        yo = matrizLocal[i][j];
        e_arriba = matrizLocal[i - 1][j];
        e_abajo = fila_abajo[j];
        e_izq = matrizLocal[i][j - 1];
        e_der = columna_derecha[i];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);

        // Se intercambian los punteros para que matrizLocal contenga los valores de temperatura recientemente calculados
        aux = matrizLocal;
        matrizLocal = matrizSiguiente;
        matrizSiguiente = aux;

    } // Fin de pasos

    time_spent = sampleTime() - time_spent;
    printf("Tlado: %d, Pasos: %d\n", Tlado, pasos);
    printf("Tiempo de ejecución: %.5f segundos.\n", time_spent);

    // Se almacena la matriz final en un archivo
    char nombre[30];
    i = j = 0;
    sprintf(nombre, "subgrid_%d_%d.out", i, j);
    FILE *f = fopen(nombre, "w");
    if (f == NULL)
    {
        printf("ERROR: No se pudo abrir el archivo\n");
        exit(1);
    }
    for (i = 0; i < Tlado; i++)
    {
        for (j = 0; j < Tlado; j++)
            fprintf(f, "%8.3f ", matrizLocal[i][j]);
        fprintf(f, "\n");
    }
    fclose(f);

    return 0;
}
