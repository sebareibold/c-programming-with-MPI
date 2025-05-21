#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>

#define Cx 0.1f
#define Cy 0.1f

double sampleTime() {
   struct timespec tv;
   clock_gettime(CLOCK_MONOTONIC,&tv);
   return ((double)tv.tv_sec+((double)tv.tv_nsec)/1000000000.0);
}

int main(int argc, char *argv[]) {
   register int p, i, j;
   if(argc != 3) {
      printf("Parámetros: Tlado pasos\n");
      exit(1);
   }
   int Tlado = atoi(argv[1]); //N --> Tamaño de la matriz cuadrada
   int pasos = atoi(argv[2]); //M --> Cantidad de pasos
   if (Tlado <= 0 || pasos < 0) {
      printf("ERROR: valores incorrectos en parámetros\n");
      exit(1);
   }

   int columnas = Tlado;
   int filas = Tlado;
   int *fila_arriba = malloc(sizeof(float)*columnas);
   int *fila_abajo = malloc(sizeof(float)*columnas);
   int *columna_izquierda = malloc(sizeof(float)*filas);
   int *columna_derecha = malloc(sizeof(float)*filas);
   bzero(fila_arriba, sizeof(float)*columnas);
   bzero(fila_abajo, sizeof(float)*columnas);
   bzero(columna_izquierda, sizeof(float)*filas);
   bzero(columna_derecha, sizeof(float)*filas);

   float **matrizActual;
   float **matrizSiguiente;
   float **aux; //Puntero auxiliar para intercambio de punteros

   //Reserva de espacio para la matriz Actual
   matrizActual = (float **)malloc(sizeof(float *) * Tlado);
   *matrizActual = (float *)malloc(sizeof(float) * Tlado * Tlado);

   //Reserva de espacio para la matriz Siguiente
   matrizSiguiente  = (float **)malloc(sizeof(float *) * Tlado);
   *matrizSiguiente = (float *)malloc(sizeof(float) * Tlado * Tlado);

   if (matrizActual == NULL || matrizSiguiente == NULL){
      printf("ERROR: No se pudo reservar memoria\n");
      exit(2);
   }

   //Se calculan las direcciones de las filas
   for(i = 1; i < Tlado; i++){
      matrizActual[i] = *matrizActual + Tlado * i;
      matrizSiguiente[i] = *matrizSiguiente + Tlado * i;
   }

   //Inicialización de la matriz Actual
   for (i = 0; i < Tlado; i++)
      for (j = 0; j < Tlado; j++)
         matrizActual[i][j] = (float)((i+1)) * (Tlado + i) * (j+1) * (Tlado + j);

   //Ejecución de los pasos de simulación
   double time_spent = sampleTime();

   float e_arriba, e_abajo, e_izq, e_der, yo;
   for (p = 0; p < pasos; p++) {

      //Se procesa el interior
      for (i = 1; i < filas-1; i++)
         for (j = 1; j < columnas-1; j++) {
            yo = matrizActual[i][j];
            e_arriba = matrizActual[i-1][j];
            e_abajo = matrizActual[i+1][j];
            e_izq = matrizActual[i][j-1];
            e_der = matrizActual[i][j+1];
            matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);
         }
      
      //Se procesa fila superior
      i = 0;
      for (j = 1; j < columnas-1; j++) {
         yo = matrizActual[i][j];
         e_arriba = fila_arriba[j];
         e_abajo = matrizActual[i+1][j];
         e_izq = matrizActual[i][j-1];
         e_der = matrizActual[i][j+1];
         matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);

      }
      
      //Se procesa fila inferior
      i = filas - 1;
      for (j = 1; j < columnas-1; j++) {
         yo = matrizActual[i][j];
         e_arriba = matrizActual[i-1][j];
         e_abajo = fila_abajo[j];
         e_izq = matrizActual[i][j-1];
         e_der = matrizActual[i][j+1];
         matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);
       }
      
      //Se procesa columna izquierda
      j = 0;
      for (i = 1; i < filas-1; i++) {
         yo = matrizActual[i][j];
         e_arriba = matrizActual[i-1][j];
         e_abajo = matrizActual[i+1][j];
         e_izq = columna_izquierda[i];
         e_der = matrizActual[i][j+1];
         matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);
      }
      
      //Se procesa columna derecha
      j = columnas - 1;
      for (i = 1; i < filas-1; i++) {
         yo = matrizActual[i][j];
         e_arriba = matrizActual[i-1][j];
         e_abajo = matrizActual[i+1][j];
         e_izq = matrizActual[i][j-1];;
         e_der = columna_derecha[i];
         matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);
      }

      //Se procesa esquina superior izquierda
      i = 0; j = 0;
      yo = matrizActual[i][j];
      e_arriba = fila_arriba[0];
      e_abajo = matrizActual[i+1][j];
      e_izq = columna_izquierda[0];
      e_der = matrizActual[i][j+1];
      matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);
       
      //Se procesa esquina superior derecha
      i = 0; j = columnas-1;
      yo = matrizActual[i][j];
      e_arriba = fila_arriba[columnas-1];
      e_abajo = matrizActual[i+1][j];
      e_izq = matrizActual[i][j-1];
      e_der = columna_derecha[i];
      matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);

      //Se procesa esquina inferior izquierda
      i = filas-1; j = 0;
      yo = matrizActual[i][j];
      e_arriba = matrizActual[i-1][j];
      e_abajo = fila_abajo[j];
      e_izq = columna_izquierda[filas-1];
      e_der = matrizActual[i][j+1];
      matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);

      //Se procesa esquina inferior derecha
      i = filas-1; j = columnas-1;
      yo = matrizActual[i][j];
      e_arriba = matrizActual[i-1][j];
      e_abajo = fila_abajo[j];
      e_izq = matrizActual[i][j-1];
      e_der = columna_derecha[i];
      matrizSiguiente[i][j] = yo+Cx*(e_abajo+e_arriba-2*yo)+Cy*(e_der+e_izq-2*yo);

      //Se intercambian los punteros para que matrizActual contenga los valores de temperatura recientemente calculados
      aux = matrizActual;
      matrizActual = matrizSiguiente;
      matrizSiguiente = aux;
   } //Fin de pasos

   time_spent = sampleTime() - time_spent;
   printf("Tlado: %d, Pasos: %d\n", Tlado, pasos);
   printf("Tiempo de ejecución: %.5f segundos.\n", time_spent);

   //Se almacena la matriz final en un archivo
   char nombre[30];
   i = j = 0;
   sprintf(nombre, "subgrid_%d_%d.out", i, j);
   FILE *f = fopen(nombre, "w");
   if (f == NULL) {
      printf("ERROR: No se pudo abrir el archivo\n");
      exit(1);
   }
   for (i = 0; i <  Tlado; i++) {
      for (j = 0; j < Tlado; j++)
         fprintf(f, "%8.3f ", matrizActual[i][j]);
      fprintf(f, "\n");
   }
   fclose(f);

   return 0;
}

