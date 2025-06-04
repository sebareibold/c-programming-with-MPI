#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>
#include <mpi.h>
#include <math.h>

#define Cx 0.1f
#define Cy 0.1f

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

    /* ======================================================= CREACION DEL MAPA: ============================================================================ */

    MPI_Comm COMM_CART;
    MPI_Status status;

    int rank, // RANK del proceso en el COMM_WORLD
        size, // Cantidad de procesos en el COMM_WORLD
        rank_cart, // RANK del proceso en el COMM_CART
        arriba, // RANK del vecino de arriba
        abajo,  // RANK del vecino de abajo
        izq,    // RANK del vecino de la izquierda
        der;    // RANK del vecino de la derecha

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Obtenemos el rank del proceso en el COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Obtenemos la cantidad de procesos en el COMM_WORLD

    int dims[2]; // Dimensiones del mapa (2D) 
    calcularDIms(size, dims); // Calculamos las dimensiones de la topologia cartesiana

    printf("Mapa [%d][%d], con size: [%d] \n", dims[0], dims[1], size);

    int periods[] = {0, 0}; // No sera periodica ninguna dimension.

    int coordenadas[2]; // Coordenadas del proceso en el COMM_CART

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &COMM_CART);
    // Creamos el COMM_CART, que es un comunicador cartesiano de 2 dimensiones (no periodico)

    MPI_Comm_rank(COMM_CART, &rank_cart);  // Obtenemos el rank del proceso en el COMM_CART
    MPI_Cart_coords(COMM_CART, rank_cart, 2, coordenadas); // Obtenemos las coordenadas del proceso en el COMM_CART

    // Calculamos la cantidad de filas y columnas segun las dimensiones de la topologia cartesiana

    int filas = Tlado / dims[0]; // Division entera filas
    int columnas = Tlado / dims[1]; // Division entera columnas


    int exceso_fila = Tlado % dims[0];
    int exceso_columna = Tlado % dims[1]; 

    int local_filas = filas + (coordenadas[0] < exceso_fila ? 1 : 0);   // Si la coordenada es menor al exceso, entonces se agrega una fila mas
    int local_columnas = columnas + (coordenadas[1] < exceso_columna ? 1 : 0); // Si la coordenada es menor al exceso, entonces se agrega una columna mas

    // El offset es la cantidad de filas y columnas que se deben saltar para llegar a la posicion del proceso en la matriz global
    // Si bien cada proceso trabaja con una matriz de dimensiones locales, el offset es necesario para inicializar correctamente los valores de la matriz local
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

    // Se reservan memoria para las matrices locales y siguientes
    // matrizLocal es la matriz que contiene los valores de temperatura en el paso actual
    // matrizSiguiente es la matriz que contiene los valores de temperatura en el paso siguiente
    // aux es una variable auxiliar para intercambiar los punteros de las matrices
    float **matrizLocal;
    float **matrizSiguiente;
    float **aux;

    matrizLocal = (float **)malloc(sizeof(float *) * local_filas);
    *matrizLocal = (float *)malloc(sizeof(float) * local_columnas * local_filas);
    matrizSiguiente = (float **)malloc(sizeof(float *) * local_filas);
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
        
        bool esi = false, esd= false, eii = false, eid = false;
        
        // Obtenemos los RANKS de los vecinos de arriba y abajo
        MPI_Cart_shift(COMM_CART, 0, 1, &arriba, &abajo); 

        // lo mismo para la dimension 1,RANK de los vecinos izquierda y derecha
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
            
            if (der == MPI_PROC_NULL) {
            	i = 0;
		j = local_columnas - 1;
		yo = matrizLocal[i][j];
		e_arriba = fila_arriba[local_columnas - 1];
		e_abajo = matrizLocal[i + 1][j];
		e_izq = matrizLocal[i][j - 1];
		e_der = columna_derecha[i];
		matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
		esd = true;
	   }
	   
	   if (izq == MPI_PROC_NULL) {
	        i = 0;
		j = 0;
		yo = matrizLocal[i][j];
		e_arriba = fila_arriba[0];
		e_abajo = matrizLocal[i + 1][j];
		e_izq = columna_izquierda[0];
		e_der = matrizLocal[i][j + 1];
		matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
		esi = true;
	   
	   }
	   bool franja_sup_procesada = true;
           
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
            
            if (der == MPI_PROC_NULL) {
            	i = local_filas - 1;
		j = local_columnas - 1;
		yo = matrizLocal[i][j];
		e_arriba = matrizLocal[i - 1][j];
		e_abajo = fila_abajo[j];
		e_izq = matrizLocal[i][j - 1];
		e_der = columna_derecha[i];
		matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
		eid = true;
	   }
	   
	   if (izq == MPI_PROC_NULL) {
	        i = local_filas - 1;
		j = 0;
		yo = matrizLocal[i][j];
		e_arriba = matrizLocal[i - 1][j];
		e_abajo = fila_abajo[j];
		e_izq = columna_izquierda[local_filas - 1];
		e_der = matrizLocal[i][j + 1];
		matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
		eii = true;	   
	   }
	   bool franja_inf_procesada = true;
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
            bool franja_izq_procesada = true;
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
            bool franja_der_procesada = true;
        }

        int recepciones_completadas = 0; // Contador de recepciones completadas
        
        while (recepciones_completadas < count_real) 
        {
            
            int indice_completado; 
            
            // Esperamos a que se complete alguna de las recepciones 
            MPI_Waitany(4, array_request_recv, &indice_completado, &status_recv_waitany); 
            recepciones_completadas++;
            
            // Procesamos la recepción completada
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
                    if (!eii && !franja_izq_procesada) {
                    	i = local_filas - 1;
			j = 0;
			yo = matrizLocal[i][j];
			e_arriba = matrizLocal[i - 1][j];
			e_abajo = fila_abajo[j];
			e_izq = columna_izquierda[local_filas - 1];
			e_der = matrizLocal[i][j + 1];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    if (!eid && !franja_der_procesada) {
                    	i = local_filas - 1;
			j = local_columnas - 1;
			yo = matrizLocal[i][j];
			e_arriba = matrizLocal[i - 1][j];
			e_abajo = fila_abajo[j];
			e_izq = matrizLocal[i][j - 1];
			e_der = columna_derecha[i];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
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
                    if (!esi && !franja_izq_procesada) {
                    	i = 0;
			j = 0;
			yo = matrizLocal[i][j];
			e_arriba = fila_arriba[0];
			e_abajo = matrizLocal[i + 1][j];
			e_izq = columna_izquierda[0];
			e_der = matrizLocal[i][j + 1];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    if (!esd && !franja_der_procesada) {
                    	i = 0;
			j = local_columnas - 1;
			yo = matrizLocal[i][j];
			e_arriba = fila_arriba[local_columnas - 1];
			e_abajo = matrizLocal[i + 1][j];
			e_izq = matrizLocal[i][j - 1];
			e_der = columna_derecha[i];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
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
                    if (!eid && !franja_inf_procesada) {
                    	i = local_filas - 1;
			j = local_columnas - 1;
			yo = matrizLocal[i][j];
			e_arriba = matrizLocal[i - 1][j];
			e_abajo = fila_abajo[j];
			e_izq = matrizLocal[i][j - 1];
			e_der = columna_derecha[i];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                    if (!esd && !franja_sup_procesada) {
                    	i = 0;
			j = local_columnas - 1;
			yo = matrizLocal[i][j];
			e_arriba = fila_arriba[local_columnas - 1];
			e_abajo = matrizLocal[i + 1][j];
			e_izq = matrizLocal[i][j - 1];
			e_der = columna_derecha[i];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
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
                    if (!eii && !franja_inf_procesada) {
                    	i = local_filas - 1;
			j = 0;
			yo = matrizLocal[i][j];
			e_arriba = matrizLocal[i - 1][j];
			e_abajo = fila_abajo[j];
			e_izq = columna_izquierda[local_filas - 1];
			e_der = matrizLocal[i][j + 1];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
				    }
                    if (!esi && !franja_sup_procesada) {
                    	i = 0;
			j = 0;
			yo = matrizLocal[i][j];
			e_arriba = fila_arriba[0];
			e_abajo = matrizLocal[i + 1][j];
			e_izq = columna_izquierda[0];
			e_der = matrizLocal[i][j + 1];
			matrizSiguiente[i][j] = yo + Cx * (e_abajo + e_arriba - 2 * yo) + Cy * (e_der + e_izq - 2 * yo);
                    }
                }
                break;
            default:
                break;
            }
        }


        MPI_Waitall(4, array_request_send, statuses_send_waitall); 

        // Se intercambian los punteros para que matrizLocal contenga los valores de temperatura recientemente calculados
        aux = matrizLocal;
        matrizLocal = matrizSiguiente;
        matrizSiguiente = aux;

    }

    // time_spent = sampleTime() - time_spent;
    printf("Tlado: %d, Pasos: %d\n", Tlado, pasos);
    // printf("Tiempo de ejecución: %.5f segundos.\n", time_spent);

    // Se almacena la matriz final en un archivo
    char nombre[30];
    i = j = 0;
    int coords[2];
    MPI_Cart_coords(COMM_CART, rank_cart, 2, coords);

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
