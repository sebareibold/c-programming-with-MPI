#set xrange [0:200] #Rango del eje X
#set yrange [10:135] #Rango del eje Y
#set xtic 5 #Valor del paso o tic del eje X
#set ytic 5 #Valor de paso o tic del eje Y

set title "Heat Map"
set view map
set xlabel "Nombre eje X"
set ylabel "Nombre eje Y"
set term png
#Válido a partir de la versión 5. En versiones previas colocar los nombres en el script.
datos=ARG1;
salida=ARG2;

#Guarda el gráfico en una imagen
set output salida.'.png'
splot datos matrix with image
