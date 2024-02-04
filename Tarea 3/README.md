# Tarea 3 Métodos Numéricos 
## Nombre: Rafael Alejandro García Ramírez 
## email: rafael.ramirez@cimat.mx
### Descripción general 
En el siguiente documento se presentan las características y funcionamiento de los códigos realizados durante la tarea 3 así como el procedimiento para su utilización.

### Uso del código : Compilación y ejecución

Los siguientes códigos fueron creadosy probados utilizando el lenguaje de programación C, más en específico con la versión **_gcc (GCC) 13.2.1 20230728 (Red Hat 13.2.1-1)_**. El nombre del archivo para cada problema consta del formate "problema **#** .c" donde **#** corresponde al número del problema que se desea ejecutar. Un ejemplo es el siguiente utilizando **bash**, que deberá ejecutarse una vez estando dentro del directorio del archivo:

`gcc problema1.c -o problema1 && ./problema1`

Donde la parte antes del `&&` corresponde a la compilación y la segunda parte a la ejecución.

## Problema 1
El perteneciente al archivo `programa1.c`. Solo se requiere del uso de la paquetería `<stdio.h>`.
### Documentación
El código proporcionado es un programa en C que resuelve un sistema de ecuaciones lineales de la forma $Lx = b$, donde $L$ es una matriz triangular inferior, $x$ es el array de incógnitas que se desea calcular y $b$ es el array de términos independientes.

#### Variables y Funciones

A continuación, se detallan las variables y funciones utilizadas en el código:

- `int N`: Esta variable representa el tamaño de la matriz $L$ y los arrays $b$ y $x$. En el ejemplo proporcionado, $N$ está definido como 3, lo que significa que estamos resolviendo un sistema de ecuaciones con 3 incógnitas.

- `double L[N][N]`: Esta es una matriz de tamaño $N\times N$ que representa la matriz triangular inferior $L$ del sistema de ecuaciones. En el ejemplo, se inicializa con valores específicos al problema dado.

- `double b[N]`: Este es un array de tamaño $N$ que representa los términos independientes del sistema de ecuaciones. En el ejemplo, se inicializa con valores específicos al problema dado.

- `double x[N]`: Este es el array de incógnitas que se desea calcular. Se declara como un array de tamaño $N$, pero su contenido se calcula utilizando la función `Lx_b`.

- `void Lx_b(int N, double L[N][N], double b[N], double x[N])`: Esta es una función que realiza el cálculo de la solución del sistema de ecuaciones $Lx = b$ utilizando el método de sustitución hacia adelante (Forward Substitution). Toma como entrada el tamaño $N$, la matriz triangular inferior $L$, el array de términos independientes $b$ y el array de incógnitas $x$. Calcula los valores de $x$ y los almacena dentro del mismo array $x$.

#### Formato de las Variables de Entrada

El formato de las variables de entrada es el siguiente:

- `N` es un número entero que representa el tamaño de la matriz y los arrays.
- `L` es una matriz de tamaño $N \times N$ de tipo double.
- `b` es un array de tamaño $N$ de tipo double.
- `x` es un array de tamaño $N$ de tipo double.

#### Funcionamiento del Programa

El programa resuelve el sistema de ecuaciones lineales $Lx = b$ utilizando el método de sustitución hacia adelante. El proceso se realiza de la siguiente manera:

1. Se define la matriz triangular inferior $L$ y el array de términos independientes $b$.

2. Se llama a la función `Lx_b` pasando como argumentos el tamaño $N$, la matriz $L$, el array $b$ y el array $x$.

3. La función `Lx_b` calcula los valores de $x$ utilizando la sustitución hacia adelante y los almacena en el array $x$.

4. El programa imprime los valores de $x$ calculados.

#####  Modificación para Entradas Más Grandes

Si se desea resolver sistemas de ecuaciones con un tamaño mayor, se puede realizar una simple modificación dentro de la función `int main()` para la lectura de un archivo o para la introducción de los valores en un `scanf()`, así como con sus respectivas entradas. 
### Código

## Problema 2
### Documentación

### Código

## Problema 3 
### Documentación

### Código

## Problema 4 
### Documentación

### Código

