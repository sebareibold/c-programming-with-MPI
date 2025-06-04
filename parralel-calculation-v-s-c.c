#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>
#include <mpi.h>
#include <math.h>

#define Cx 0.1f
#define Cy 0.1f
/*

    1) Cada proceso hace la inicializacion de su MATRIZ (TENE encuenta las coordenadasenanas del MAPA y Tlado)
    2) Cominucacio(tanto recibir como enviar) --> Todo lo recibido se almacena en los vectores auxiliares
    3) Ya con los datos hacemos el calculo --> Imprimimos nuestro .out
    4) El script genera el .out general a partir de los otro sub.... .out
*/

double sampleTime()
{
    struct timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return ((double)tv.tv_sec + ((double)tv.tv_nsec) / 1000000000.0);
}

void calcularDIms(int size, int dims[2])
{
    dims[0] = 1;
    dims[1] = size;
    for (int i = 1; i <= sqrt(size); i++)
    {
        if (size % i == 0)
        {
            dims[0] = i;
            dims[1] = size / i;
        }
    }
}
int main(int argc, char *argv[])
{
    int p, i, j;
    if (argc != 3)
    {
        printf("Parámetros: Tlado pasos\n");
        exit(1);
    }

    MPI_Init(&argc, &argv);

    int Tlado = atoi(argv[1]); // N --> Tamaño de la matriz cuadrada
    int pasos = atoi(argv[2]); // M --> Cantidad de pasos

    if (Tlado <= 0 || pasos < 0)
    {
        printf("ERROR: valores incorrectos en parámetros\n");
        exit(1);
    }

    /* ============================================================= DECLARACION DE VARIABLES  =============================================================== */
    MPI_Comm COMM_CART;
    MPI_Status status;

    /* ======================================================= CREACION DEL MAPA: 2 DIMENSIONES (NO PERIODICO) ========================================================== */

    int rank, size, rank_cart, arriba, abajo, izq, der;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[2];
    calcularDIms(size, dims);
    printf("Mapa [%d][%d], con size: [%d] \n", dims[0], dims[1], size);

    int periods[] = {0, 0}; // No sera periodica ninguna dimension.

    int coordenadas[2];

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &COMM_CART);
    MPI_Comm_rank(COMM_CART, &rank_cart);                  // Obtenemos el rank en el nuevo Communicador
    MPI_Cart_coords(COMM_CART, rank_cart, 2, coordenadas); // Obtenemos las coordenas

    int cant_filas_mapa = dims[0], cant_columnas_mapa = dims[1];

    // Caso TLado Divisible (por cantidad de procesos)
    int filas = Tlado / dims[0];
    int columnas = Tlado / dims[1];
    printf("Soy Proceso %d \n", rank_cart);
    printf("Filas (div de Tlado): %d \n", filas);
    printf("Columnas (div de Tlado): %d\n", columnas);

    // Tenes encuenta en resto de la division de arriba
    int exceso_fila = Tlado % dims[0];
    int exceso_columna = Tlado % dims[1];

    // Obtenemos las filas locales de la matriz, teniendo en cuenta si hubo un exceso, si lo hubo entonces le agregamos una fila o columnas (por las 2 lineas)
    int local_filas = filas + (coordenadas[0] < exceso_fila ? 1 : 0);
    int local_columnas = columnas + (coordenadas[1] < exceso_columna ? 1 : 0);
    printf("Local filas: %d\n", local_filas);
    printf("Local Columnas: %d\n", local_columnas);

    // Inicializa su inicio de fila y columna a partir de las coordenadas y fijas (o columnas), para despues si es que el exceso es menor a la coordenada entonces coordenas
    int offset_fila = coordenadas[0] * filas + (coordenadas[0] < exceso_fila ? coordenadas[0] : exceso_fila);
    int offset_columna = coordenadas[1] * columnas + (coordenadas[1] < exceso_columna ? coordenadas[1] : exceso_columna);

    /* =======================================================  SECCION RELACIONADA A MEMORIA (MATRIZ Y ARREGLOS) ======================================================= */

    // Bloque addiciones para cada el envio de filas o columnas vecinas
    float *fila_arriba = (float *)malloc(sizeof(float) * local_columnas);
    float *fila_abajo = (float *)malloc(sizeof(float) * local_columnas);
    float *columna_izquierda = (float *)malloc(sizeof(float) * local_filas);
    float *columna_derecha = (float *)malloc(sizeof(float) * local_filas);

    bzero(fila_arriba, sizeof(float) * local_columnas);
    bzero(fila_abajo, sizeof(float) * local_columnas);
    bzero(columna_izquierda, sizeof(float) * local_filas);
    bzero(columna_derecha, sizeof(float) * local_filas);

    // Matriz Local de cada Proceso, dicha matriz es de TLado X TLado
    float **matrizLocal;
    float **matrizSiguiente;
    float **aux;

    matrizLocal = (float **)malloc(sizeof(float *) * local_filas); // Declaracion eficiente de la matriz
    *matrizLocal = (float *)malloc(sizeof(float) * local_columnas * local_filas);
    matrizSiguiente = (float **)malloc(sizeof(float *) * local_filas); // Declaracion eficiente de la matriz
    *matrizSiguiente = (float *)malloc(sizeof(float) * local_filas * local_columnas);

    if (matrizLocal == NULL || matrizSiguiente == NULL)
    {
        printf("ERROR: No se pudo reservar memoria\n");
        exit(2);
    }

    for (i = 1; i < local_filas; i++)
    {
        matrizLocal[i] = *matrizLocal + local_columnas * i;
        matrizSiguiente[i] = *matrizSiguiente + local_columnas * i;
    }

    /* ============================================================= (1) INICIALIZACION DE LA MATRIZ  =============================================================== */

    for (i = 0; i < local_filas; i++)
        for (j = 0; j < local_columnas; j++)
        {
            matrizLocal[i][j] = (float)((i + offset_fila + 1)) * (Tlado + i + offset_fila) * (j + offset_columna + 1) * (Tlado + j + offset_columna);
        }

    /* ======================================================= (2) SECCION DE ENVIO Y RECEPCION CON LOS VECINOS ========================================================== */

    // Creamos el tipo de datos para pasar las columnas de la derecha e izquierda a los vecinos
    MPI_Datatype vectorVertical;
    MPI_Type_vector(local_filas, 1, local_columnas, MPI_FLOAT, &vectorVertical);
    MPI_Type_commit(&vectorVertical);

    
    MPI_Request array_request_recv[4]; 
    MPI_Request array_request_send[4];
     
    MPI_Status status_recv_waitany; 
    MPI_Status statuses_send_waitall[4];
    
    int count = 4, count_real;
    
    float e_arriba, e_abajo, e_izq, e_der, yo;
    for (p = 0; p < pasos; p++)
    {
        // Se inicializan los requests a MPI_REQUEST_NULL en cada paso.
        for(int k=0; k<4; ++k) {
            array_request_recv[k] = MPI_REQUEST_NULL;
            array_request_send[k] = MPI_REQUEST_NULL;
        }
 
        count_real = 0;
        
        // Variables auxiliares para procesar las franjas (para no pasar por alto los bordes fisicos).
        bool franja_sup_procesada = false;
        bool franja_inf_procesada = false;
        bool franja_izq_procesada = false;
        bool franja_der_procesada = false;
        
        // Obtenemos los RANKS de los vecinos de arriba y abajo
        MPI_Cart_shift(COMM_CART, 0, 1, &arriba, &abajo); // dame los vecinos (destino)

        // lo mismo para la dimension 1,RANK de los vecinos izq y der
        MPI_Cart_shift(COMM_CART, 1, 1, &izq, &der);

        if (arriba != MPI_PROC_NULL) // si tengo vecino arriba MANDO mi fila superior
        {
            MPI_Isend(&matrizLocal[0][0], local_columnas, MPI_FLOAT, arriba, 0, COMM_CART, &array_request_send[0]);
        }
        if (abajo != MPI_PROC_NULL) // si tengo vecino abajo RECIBO su fila superior
        {
            MPI_Irecv(&fila_abajo[0], local_columnas, MPI_FLOAT, abajo, 0, COMM_CART, &array_request_recv[0]);
            count_real++;
        }
        if (abajo != MPI_PROC_NULL) // si tengo vecino abajo MANDO mi fila inferior
        {
            MPI_Isend(&matrizLocal[local_filas - 1][0], local_columnas, MPI_FLOAT, abajo, 0, COMM_CART, &array_request_send[1]);
        }
        if (arriba != MPI_PROC_NULL) // si tengo vecino arriba RECIBO su fila inferior
        {
            MPI_Irecv(&fila_arriba[0], local_columnas, MPI_FLOAT, arriba, 0, COMM_CART, &array_request_recv[1]);
            count_real++;
        }
        if (izq != MPI_PROC_NULL) // si tengo vecino izquierda mando mi columna IZQ.
        {
            MPI_Isend(&matrizLocal[0][0], 1, vectorVertical, izq, 0, COMM_CART, &array_request_send[2]);
        }
        if (der != MPI_PROC_NULL) // si tengo vecino derecha RECIBO su columna IZQ (seria MI COLUMNA EXTERIOR UBICADA A LA DERECHA).
        {
            MPI_Irecv(&columna_derecha[0], local_filas, MPI_FLOAT, der, 0, COMM_CART, &array_request_recv[2]);
            count_real++;
        }

        if (der != MPI_PROC_NULL) // si tengo vecino derecha mando mi columna derecha
        {
            MPI_Isend(&matrizLocal[0][local_columnas - 1], 1, vectorVertical, der, 0, COMM_CART, &array_request_send[3]);
        }
        if (izq != MPI_PROC_NULL) // si tengo vecino izquierda RECIBO su columna derecha (seria MI COLUMNA EXTERIOR UBICADA A LA IZQUIERDA).
        {
            MPI_Irecv(&columna_izquierda[0], local_filas, MPI_FLOAT, izq, 0, COMM_CART, &array_request_recv[3]);
            count_real++;
        }

        /* ======================================================= (3) CALCULO DE LAS SECCION INTERIOR Y DE LOS BORDES ========================================================== */

        // Se procesa el interior
        for (i = 1; i < local_filas - 1; i++)
            for (j = 1; j < local_columnas - 1; j++)
            {
                yo = matrizLocal[i][j];
                e_arriba = matrizLocal[i - 1][j];
                e_abajo = matrizLocal[i + 1][j];
                e_izq = matrizLocal[i][j - 1];
                e_der = matrizLocal[i][j + 1];
                matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
            }

        int recepciones_completadas = 0;    
        
        while (recepciones_completadas < count_real) 
        {
            int indice_completado; 
            MPI_Waitany(4, array_request_recv, &indice_completado, &status_recv_waitany); 
            recepciones_completadas++;
            
            switch (indice_completado)
            {
            case 0:
                if (!franja_inf_procesada && local_filas > 1) { 
                    i = local_filas - 1;
                    for (j = 1; j < local_columnas - 1; j++) {
                        yo = matrizLocal[i][j];
                        e_arriba = matrizLocal[i - 1][j];
                        e_abajo = fila_abajo[j];
                        e_izq = matrizLocal[i][j - 1];
                        e_der = matrizLocal[i][j + 1];
                        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    franja_inf_procesada = true;
                }
                break;
            case 1:
                if (!franja_sup_procesada && local_filas > 0) {
                    i = 0;
                    for (j = 1; j < local_columnas - 1; j++) {
                        yo = matrizLocal[i][j];
                        e_arriba = fila_arriba[j];
                        e_abajo = (local_filas > 1) ? matrizLocal[i + 1][j] : fila_abajo[j];
                        e_izq = matrizLocal[i][j - 1];
                        e_der = matrizLocal[i][j + 1];
                        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    franja_sup_procesada = true;
                }
                break;
            case 2:
                if (!franja_der_procesada && local_columnas > 1) {
                    j = local_columnas - 1;
                    for (i = 1; i < local_filas - 1; i++) {
                        yo = matrizLocal[i][j];
                        e_arriba = matrizLocal[i - 1][j];
                        e_abajo = matrizLocal[i + 1][j];
                        e_izq = matrizLocal[i][j - 1];
                        e_der = columna_derecha[i];
                        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    franja_der_procesada = true;
                }
                break;
            case 3: 
                if (!franja_izq_procesada && local_columnas > 0) {
                    j = 0;
                    for (i = 1; i < local_filas - 1; i++) {
                        yo = matrizLocal[i][j];
                        e_arriba = matrizLocal[i - 1][j];
                        e_abajo = matrizLocal[i + 1][j];
                        e_izq = columna_izquierda[i];
                        e_der = (local_columnas > 1) ? matrizLocal[i][j + 1] : columna_derecha[i];
                        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    franja_izq_procesada = true;
                }
                break;
            default:
                break;
            }
        }

        if (arriba == MPI_PROC_NULL && !franja_sup_procesada && local_filas > 0) {
            i = 0;
            for (j = 1; j < local_columnas - 1; j++) { 
                yo = matrizLocal[i][j];
                e_arriba = fila_arriba[j];
                e_abajo = (local_filas > 1) ? matrizLocal[i + 1][j] : fila_abajo[j];
                e_izq = matrizLocal[i][j - 1];
                e_der = matrizLocal[i][j + 1];
                matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2*yo) + Cy * (e_der + e_izq - 2*yo);
            }
           
        }
        
        if (abajo == MPI_PROC_NULL && !franja_inf_procesada && local_filas > 1) {
            i = local_filas - 1;
            for (j = 1; j < local_columnas - 1; j++) {
                yo = matrizLocal[i][j];
                e_arriba = matrizLocal[i - 1][j];
                e_abajo = fila_abajo[j];
                e_izq = matrizLocal[i][j - 1];
                e_der = matrizLocal[i][j + 1];
                matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2*yo) + Cy * (e_der + e_izq - 2*yo);
            }
        }
        if (izq == MPI_PROC_NULL && !franja_izq_procesada && local_columnas > 0) {
            j = 0;
            for (i = 1; i < local_filas - 1; i++) {
                yo = matrizLocal[i][j];
                e_arriba = matrizLocal[i - 1][j];
                e_abajo = matrizLocal[i + 1][j];
                e_izq = columna_izquierda[i]; 
                e_der = (local_columnas > 1) ? matrizLocal[i][j + 1] : columna_derecha[i]; 
                matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2*yo) + Cy * (e_der + e_izq - 2*yo);
            }
        }
        if (der == MPI_PROC_NULL && !franja_der_procesada && local_columnas > 1) {
            j = local_columnas - 1;
            for (i = 1; i < local_filas - 1; i++) {
                yo = matrizLocal[i][j];
                e_arriba = matrizLocal[i - 1][j];
                e_abajo = matrizLocal[i + 1][j];
                e_izq = matrizLocal[i][j - 1];
                e_der = columna_derecha[i];
                matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2*yo) + Cy * (e_der + e_izq - 2*yo);
            }
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
        j = local_columnas - 1;
        yo = matrizLocal[i][j];
        e_arriba = fila_arriba[local_columnas - 1];
        e_abajo = matrizLocal[i + 1][j];
        e_izq = matrizLocal[i][j - 1];
        e_der = columna_derecha[i];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);

        // Se procesa esquina inferior izquierda
        i = local_filas - 1;
        j = 0;
        yo = matrizLocal[i][j];
        e_arriba = matrizLocal[i - 1][j];
        e_abajo = fila_abajo[j];
        e_izq = columna_izquierda[local_filas - 1];
        e_der = matrizLocal[i][j + 1];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);

        // Se procesa esquina inferior derecha
        i = local_filas - 1;
        j = local_columnas - 1;
        yo = matrizLocal[i][j];
        e_arriba = matrizLocal[i - 1][j];
        e_abajo = fila_abajo[j];
        e_izq = matrizLocal[i][j - 1];
        e_der = columna_derecha[i];
        matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);


        MPI_Waitall(4, array_request_send, statuses_send_waitall); 

        // Se intercambian los punteros para que matrizLocal contenga los valores de temperatura recientemente calculados
        aux = matrizLocal;
        matrizLocal = matrizSiguiente;
        matrizSiguiente = aux;

    } // Fin de pasos

    // time_spent = sampleTime() - time_spent;
    printf("Tlado: %d, Pasos: %d\n", Tlado, pasos);
    // printf("Tiempo de ejecución: %.5f segundos.\n", time_spent);

    // Se almacena la matriz final en un archivo
    char nombre[30];
    i = j = 0;
    int coords[2];
    MPI_Cart_coords(COMM_CART, rank_cart, 2, coords); // Obtener coordenadas

    sprintf(nombre, "subgrid_%d_%d.out", coords[0], coords[1]); // Usar coordenadas en el nombre
    FILE *f = fopen(nombre, "w");

    if (f == NULL)
    {
        printf("ERROR: No se pudo abrir el archivo\n");
        exit(1);
    }
    for (i = 0; i < local_filas; i++)
    {
        for (j = 0; j < local_columnas; j++)
            fprintf(f, "%8.3f ", matrizLocal[i][j]);
        fprintf(f, "\n");
    }
    fclose(f);

    free(fila_arriba);
    free(fila_abajo);
    free(columna_izquierda);
    free(columna_derecha);
    free(matrizLocal);

    MPI_Type_free(&vectorVertical);

    MPI_Finalize();

    return 0;
}
