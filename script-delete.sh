#!/bin/bash

# Este script elimina archivos en Ubuntu con el patrón subgrid_X_Y.out
# donde X e Y son dígitos individuales (0-9)

FILE_PATTERN="subgrid_[0-9]_[0-9].out"

echo "Buscando y eliminando archivos con el patrón: $FILE_PATTERN"

# --- Prueba ---
 echo "Los siguientes archivos *serían* eliminados (solo listando):"
 ls $FILE_PATTERN 2>/dev/null # 2>/dev/null evita que muestre errores si no encuentra nada

# --- Eliminación  ---
rm -v $FILE_PATTERN
 exit 0 # Salir después de la prueba

# -------------------
