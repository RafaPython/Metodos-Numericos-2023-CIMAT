/*
Tarea 6  Métodos Numéricos
Nombre : Rafael Alejandro García Ramírez
email : rafael.ramirez@cimat.mx

Ejecucion --> ./RUN M_d.txt b_d.txt TOL iteracionesMAX
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "time.h"
#include "solvers.h"


int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Error en la cantidad de argumentos.\n");
        printf("Uso: ./RUN M_d.txt b_d.txt TOL iteracionesMAX\n");
        return 1;
    }

    char nombre1[strlen(argv[1]) + 1], nombre2[strlen(argv[2])];

    double TOL = strtod(argv[3], NULL);
    int iteracionesMAX = strtod(argv[4], NULL);

    strcpy(nombre1, argv[1]);
    strcpy(nombre2, argv[2]);

    printf("Archivo matriz --> %s\n", nombre1);
    printf("Archivo vector --> %s\n", nombre2);
    printf("Tolerancia ------> %s\n", argv[3]);
    FILE *archivo_matriz;
    int filas_matriz, columnas_matriz;


    archivo_matriz = fopen(nombre1, "r");
    if (archivo_matriz == NULL) { printf("No se pudo abrir el archivo 1.\n"); return 1; }

    fscanf(archivo_matriz, "%d %d", &filas_matriz, &columnas_matriz);
    double **matriz = generar_matriz(filas_matriz);

    escritura_matriz(archivo_matriz, filas_matriz, columnas_matriz, matriz);

    FILE *archivo_vector;
    int filas_vector, columnas_vector;

    archivo_vector = fopen(nombre2, "r");
    if (archivo_vector == NULL) { printf("No se pudo abrir el archivo 2.\n"); return 1; }

    fscanf(archivo_vector, "%d %d", &filas_vector, &columnas_vector);

    double *b = generar_vector(filas_vector);
    escritura_vector(archivo_vector, filas_vector, b);
    
    int N = filas_matriz;

    printf("\nRESOLUCION  JACOBI y GAUSS-SIEDER\n");
    
    double *x0 = (double *) malloc(N * sizeof (double));

    // incializamos el vector con 1
    for (int i = 0; i < N; i++) x0[i] = 1.0 ;

if (N < 10){

    int opcion;
    printf("--- Menu ---\n");
    printf("1. Jacobi\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Salir\n");
    printf("Ingrese una opcion: ");
    scanf("%d", &opcion);

    switch(opcion){

        case 1:
            printf("Jacobi\n");
            printf("\nMatriz:\n");
            imprimir_matriz(N, matriz);
            printf("\nVector b:\n");
            imprimir_vector(N, b);
            printf("Vector x0:\n");
            imprimir_vector(N, x0);
            jacobi(N, TOL, matriz, b, x0, iteracionesMAX);
            printf("\nVector x:\n");
            imprimir_vector(N, x0);
            printf("Comprobacion: El numero de elementos diferentes a 0 del vector x es: %d\n", comprobar_solucion(N, matriz, x0, b, TOL));
            break;
        case 2:
            printf("Gauss-Seidel\n");
            printf("\nMatriz:\n");
            imprimir_matriz(N, matriz);
            printf("\nVector b:\n");
            imprimir_vector(N, b);
            printf("Vector x0:\n");
            imprimir_vector(N, x0);
            gauss_seidel(N, TOL, matriz, b, x0, iteracionesMAX);
            printf("\nVector x:\n");
            imprimir_vector(N, x0);
            printf("Comprobacion: El numero de elementos diferentes a 0 del vector x es: %d\n", comprobar_solucion(N, matriz, x0, b, TOL));
            break;
        case 3:
            printf("Salir\n");
            exit(1);
            break;
    }

} else {

    printf("\n\nLa dimension de la matriz es mayor a 10. Se guardarán la solución de salida en un archivo.\n\n");
    
    int opcion;
    printf("--- Menu ---\n");
    printf("1. Jacobi\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Salir\n");
    printf("Ingrese una opcion: ");
    scanf("%d", &opcion);

    switch(opcion){

        case 1:
            jacobi(N, TOL, matriz, b, x0, iteracionesMAX);
            printf("\nComprobacion: El numero de elementos diferentes a 0 del vector x es: %d\n", comprobar_solucion(N, matriz, x0, b, TOL));
            // creamos el nombre para el archivo con el número N 
            char nombre_archivo[100];
            sprintf(nombre_archivo, "metodo_jacobi_N%d.txt", N);

            // guardamos en un archivo
            FILE *archivo_salida;
            archivo_salida = fopen(nombre_archivo, "w");
            if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }

            for (int i = 0; i < N; i++) fprintf(archivo_salida, "%f\n", x0[i]);
            fclose(archivo_salida);
            printf("\n\nSe guardó la solución en el archivo %s\n", nombre_archivo);


            break;
        case 2:
            gauss_seidel(N, TOL, matriz, b, x0, iteracionesMAX);
            printf("\nComprobacion: El numero de elementos diferentes a 0 del vector x es: %d\n", comprobar_solucion(N, matriz, x0, b, TOL));
            // creamos el nombre para el archivo con el número N
            char nombre_archivo2[100];
            sprintf(nombre_archivo2, "metodo_gauss_seidel_N%d.txt", N);

            // guardamos en un archivo
            FILE *archivo_salida2;
            archivo_salida2 = fopen(nombre_archivo2, "w");
            if (archivo_salida2 == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }

            for (int i = 0; i < N; i++) fprintf(archivo_salida2, "%f\n", x0[i]);
            fclose(archivo_salida2);
            printf("\n\nSe guardó la solución en el archivo %s\n", nombre_archivo2);
            break;
        case 3:
            printf("Salir\n");
            exit(1);
            break;
    }
}

// liberamos la memoria
    liberar_matriz(N, matriz);
    liberar_vector(b);
    liberar_vector(x0);
    return 0;
}


