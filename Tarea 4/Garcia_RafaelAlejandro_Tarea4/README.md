
# Tarea 4: Factorización LU y Método de Cholesky
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en Linux (Fedora 38)

Resuelve sistemas de ecuaciones lineales utilizando diferentes métodos de factorización LU y el Método de Cholesky.

## Compilación

Se dividen los archivos en 2: **main.c** para los algorimos de factorización y **cholesky.c** para el método de Cholesky. Para compilar el programa, utiliza el siguiente comando según el que se requiera:

```bash
$ gcc main.c -o main -lm
$ gcc cholesky.c -o cholesky -lm
```

## Ejecución

Para ejecutar el programa principal, utiliza el siguiente comando:

```bash
$ ./main NombreMatriz.txt NombreVector.txt TOLERANCIA
```

Para ejecutar el programa con el Método de Cholesky, utiliza el siguiente comando:

```bash
$ ./cholesky
```

## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando diferentes métodos, que incluyen:

1. Método de Doolittle
2. Método de Crout
3. Método de Cholesky

Mientras que el programa cholesky.c permite solamente la resolución de la matriz de calor en 1D mediante Cholesky y guarda los resultados en un archivo .txt

### Uso y entrada

Ambos programas requieren del uso de la paquetería **_solvers.h_** dentro de la misma carpeta donde se encuentre el programa. Para utilizar el programa principal, debes proporcionar los siguientes argumentos (banderas) de línea de comandos:

- `NombreMatriz.txt`: Nombre del archivo que contiene la matriz del sistema.
- `NombreVector.txt`: Nombre del archivo que contiene el vector de términos independientes.
- `TOLERANCIA`: (**_Double_**) Valor de tolerancia para la comprobación de la solución.

El programa te dará a ecoger cuál de los tres métodos utilizar. Después, el programa realizará los siguientes pasos:

1. Lee los archivos de entrada.
2. Realiza la descomposición LU o Cholesky según el método elegido.
3. Resuelve el sistema de ecuaciones.
4. Imprime los archivos de solución en caso de que las dimensiones no sean mayores a 10, en caso contrario se guardarán los resultados en un archivo .txt 

Para el archivo de cholesky.c no se requieren mayores instrucciones que las especificadas en la sección **Ejecución**.

### Salida

Para un archivo de dimensiones menores a 10, el programa imprime cada paso de cálculo entre obtención de matrices y vectores hasta finalizar con el vetor final de solución, en caso contrario, el programa genera archivos con los resultados de la solución. 

Para el archivo cholesky.c siempre se guardarán los resultados dentro de un archivo .txt

En ambos casos los archivos generados se encontrarán dentro de la misma carperta donde se encuentre el programa.
