
# Tarea 14: Ecuaciones diferenciales parciales
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en macOS (Ventura 13.3)

Solucionador de la ecuación de calor utilizando el método $\theta$


El programa pide todas las variables que se necesiten para resolver el problema.
## Compilación

Para la compilación del archivo principal ```main.c``` se hace el comando 

```bash
$ gcc main.c -o main -lm
```

## Ejecución

Para ejecutar el programa principal, utiliza el siguiente comando con SIN banderas:

```bash
$ ./main 
```
NOTA: El programa a continuación te pedirá las constantes del problema así como sus límites. Esto lo hace para la condición inicial $4x - 4x^{2}$.
## Descripción

El programa principal resuelve para la condición incial $4x - 4x^{2}$ el problema de la ecuación de calor mediante el método $\theta$. El programa pedirá todas las variables necesarias para las condiciones específicas, así como el número de puntos y el número de iteraciones máximas. El programa resuelve el sistema de ecuaciones mediante el método de resolución de ecuaciones lineales Gauss-Seidel con una tolerancia de $\varepsilon = 1\times 10^{-5}$.

### Uso y entrada
El método te pedirá introducir las constantes necesarias para el sistema con características específicas que se desee resolver.

### Salida
Se guardará dentro de un archivo $\texttt{.txt}$ los datos generados por el método, se le pedirá al usuario decidir si desea guardar la matriz como una matriz o si desea el formato $\texttt{i,j,solucion[i][j]}$. Después imprimirá en pantalla todas las constantes que se utilizaron y un mensaje de éxito de guardado de datos. Finalmente si el usario eligió el segundo formato, se graficará usando $\texttt{gnuplot}$ los resultados.
