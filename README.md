# Simulación 2D de Difusión (Serial y MPI)

Este proyecto contiene implementaciones en lenguaje C para simular un proceso de difusión en 2D (similar a la transferencia de calor) sobre una placa cuadrada. Se incluyen dos versiones: una serial para ejecución en un solo procesador y otra paralela utilizando MPI (Message Passing Interface) para distribuir la carga computacional.

## Descripción

El código simula la evolución de una matriz cuadrada de temperaturas (o una magnitud similar) a lo largo de varios pasos de tiempo, aplicando un método de diferencias finitas (un stencil de 2D).

*   La versión **serial** (`serial-program/serial-calculation.c`) realiza la simulación completa en un único hilo de ejecución.
*   La versión **MPI** (`parralel-program/parralel-calculation-v-s-c.c`) distribuye la matriz entre múltiples procesos, utilizando una topología cartesiana 2D y comunicación entre vecinos para intercambiar los datos de frontera necesarios en cada paso de tiempo.

Ambas implementaciones toman como argumentos el tamaño del lado de la matriz cuadrada global (`Tlado`) y la cantidad de pasos de simulación (`pasos`).

La inicialización de la matriz sigue la fórmula `(i+1) * (Tlado + i) * (j+1) * (Tlado + j)` para la celda global en la posición `(i, j)`. Las condiciones de frontera externas de la placa se asumen como cero constante (condiciones de Dirichlet).

## Estructura del Proyecto

```
c-programming-with-MPI/
├── parralel-program/
│   ├── graficar_salida.gp            # Script Gnuplot para visualización
│   ├── juntar.sh                     # Script para juntar las submatrices de salida
│   ├── parralel-calculation          # Binario ejecutable (tras compilación MPI)
│   ├── parralel-calculation-v-s-c.c  # Código fuente MPI
│   ├── script-delete.sh              # Script para limpiar archivos de salida
│   └── script.sh                     # Script de utilidad (puede variar su función)
└── serial-program/
    ├── serial-calculation            # Binario ejecutable (tras compilación serial)
    └── serial-calculation.c          # Código fuente serial
```

## Requisitos

*   Un compilador de C (por ejemplo, GCC).
*   Para la versión MPI: Una implementación de MPI instalada (por ejemplo, OpenMPI, MPICH).

## Compilación

### Versión Serial

Navega al directorio `serial-program` y compila el código fuente:

```bash
cd serial-program/
gcc serial-calculation.c -o serial-calculation -lm -lrt
cd ..
```

`-lm` enlaza la librería matemática (para `sqrt`) y `-lrt` para `clock_gettime`. Aunque el tiempo de ejecución esté comentado en uno de los códigos, es buena práctica incluirlas si las funciones están presentes.

### Versión Paralela (MPI)

Navega al directorio `parralel-program` y compila el código fuente utilizando el compilador MPI:

```bash
cd parralel-program/
mpicc parralel-calculation-v-s-c.c -o parralel-calculation -lm -lrt
cd ..
```

## Ejecución

### Versión Serial

Ejecuta el binario compilado desde el directorio `serial-program`:

```bash
./serial-program/serial-calculation <Tlado> <pasos>
```

- `<Tlado>`: Tamaño del lado de la matriz cuadrada global (N).
- `<pasos>`: Cantidad de pasos de simulación (M).

Ejemplo:

```bash
./serial-program/serial-calculation 512 100
```

La matriz final se guardará en el archivo `subgrid_0_0.out` dentro del directorio donde se ejecutó el comando (usualmente `serial-program/`).

### Versión Paralela (MPI)

Ejecuta el binario compilado utilizando `mpirun` desde el directorio `parralel-program`:

```bash
mpirun -np <num_procesos> ./parralel-program/parralel-calculation <Tlado> <pasos>
```

- `<num_procesos>`: Número de procesos MPI a utilizar.
- `<Tlado>`: Tamaño del lado de la matriz cuadrada global (N).
- `<pasos>`: Cantidad de pasos de simulación (M).

Ejemplo (con 4 procesos):

```bash
mpirun -np 4 ./parralel-program/parralel-calculation 1024 200
```

Cada proceso MPI guardará su submatriz local final en un archivo separado llamado `subgrid_row_col.out` dentro del directorio donde se ejecutó el comando (usualmente `parralel-program/`), donde `row` y `col` son las coordenadas del proceso en la topología cartesiana 2D.

## Scripts Adicionales (en `parralel-program/`)

- `juntar.sh`: Script diseñado para combinar los archivos de salida `subgrid_*.out` generados por la ejecución paralela con MPI en un único archivo que represente la matriz global completa.
- `script-delete.sh`: Script para facilitar la limpieza, eliminando los archivos de salida (`subgrid_*.out`) generados por las ejecuciones.
- `graficar_salida.gp`: Un script para Gnuplot, presumiblemente configurado para visualizar los datos contenidos en los archivos de salida.
- `script.sh`: Script de utilidad general. Su función específica dependerá de su contenido, pero podría usarse para automatizar tareas de compilación, ejecución de pruebas, o post-procesamiento.

## Notas sobre el Código

- La función `calcularDIms` (utilizada en la versión MPI) tiene como objetivo encontrar las dimensiones `[dims[0]][dims[1]]` para la topología cartesiana 2D que mejor se ajusten al número de procesos (`size`), buscando una configuración lo más cercana posible a un cuadrado para optimizar la distribución de la carga y la comunicación.
- La versión MPI implementa una estrategia de descomposición de dominio en bloques 2D, donde cada proceso se encarga de simular una porción rectangular de la matriz global. El código maneja la distribución del trabajo incluso cuando el tamaño total (`Tlado`) no es perfectamente divisible por las dimensiones de la topología.
- La comunicación entre procesos vecinos para el intercambio de las filas y columnas frontera se realiza utilizando operaciones de comunicación no bloqueantes (`MPI_Isend`, `MPI_Irecv`). El uso de `MPI_Waitany` permite una superposición parcial de cálculo y comunicación.
- Las condiciones de frontera de la placa completa están fijadas a un valor de 0 en ambas implementaciones (condiciones de Dirichlet). Esto se gestiona utilizando arreglos auxiliares o valores implícitos de 0 fuera de los límites de la matriz.
