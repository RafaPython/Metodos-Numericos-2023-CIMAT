/*
Tarea 4: Métodos de descomposición LU
Nombre: Rafael Alejandro García Ramírez 
email: rafael.ramirez@cimat.mx
Ejecucion : gcc main.c -o main -lm && ./main archivoMatriz.txt archivoVector.txt TOLERANCIA
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "time.h"
#include "solvers.h"


int main(int argc, char *argv[]) {
    
    char nombre1[strlen(argv[1])], nombre2[strlen(argv[2])];

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

    double b[filas_vector];
    escritura_vector(archivo_vector, filas_vector, b);
    int N = filas_matriz;
    double TOL = strtod(argv[3], NULL);


    if (N > 10) {
        printf("La matriz es muy grande para imprimirse. El resultado será guardado en formato de archivo .txt \n");
        
        int opcion;
        printf("\nElija el método que desee aplicar:\n");    
        printf("1. Doolittle\n");
        printf("2. Crout\n");
        printf("3. Cholesky\n");
        printf("Introduzca el número de la opción: ");
        scanf("%d", &opcion);

        double **L = generar_matriz(N);
        double **U = generar_matriz(N);
        double **Lt = generar_matriz(N);
        double *x = generar_vector(N);
        double *z = generar_vector(N);
        double *Ax = generar_vector(N);

        switch (opcion) {
        case 1:
            printf("\n--- METODO DOOLITTLE ---\n");
            doolittle(N, matriz, L, U); // Descomposición LU
            Lx_b(N, L, b, z); // Resolución de Lz = b
            Ux_b(N, U, z, x); // Resolución de Ux = z
            // Comprobación Ax - b
            multiplicar_matriz_vector(N, matriz, x, Ax);
            for (int i = 0; i < N; i++) Ax[i] -= b[i];
            for (int i = 0; i < N; i++) if (fabs(Ax[i]) >= TOL) printf("La solución falló en x[%d] = %f \n", i, x[i]) ;
            
            // creamos un archivo para guardar la solución 
            FILE *archivo_solucion;
            archivo_solucion = fopen("solucion_doolittle.txt", "w");
            if (archivo_solucion == NULL) { printf("No se pudo abrir el archivo 3.\n"); return 1; }
            // escribimos la solución en el archivo
            for (int i = 0; i < N; i++) fprintf(archivo_solucion, "%f\n", x[i]); 
            fclose(archivo_solucion);
            printf("\n\nLa solución fue guardada en el archivo solucion_doolittle.txt\n");
            break;

        case 2:
            printf("\n--- METODO CROUT ---\n");
            crout(N, matriz, L, U); // Descomposición LU
            Lx_b(N, L, b, z); // Resolución de Lz = b
            Ux_b(N, U, z, x); // Resolución de Ux = z
            // Comprobación Ax - b
            multiplicar_matriz_vector(N, matriz, x, Ax);
            for (int i = 0; i < N; i++) {Ax[i] -= b[i]; if (fabs(Ax[i]) >= TOL) printf("La solución falló en x[%d] = %f \n",i, Ax[i]);}

            // creamos un archivo para guardar la solución
            archivo_solucion = fopen("solucion_crout.txt", "w");
            if (archivo_solucion == NULL) { printf("No se pudo abrir el archivo 3.\n"); return 1; }
            // escribimos la solución en el archivo
            for (int i = 0; i < N; i++) fprintf(archivo_solucion, "%f\n", x[i]);
            fclose(archivo_solucion);
            printf("\n\nLa solución fue guardada en el archivo solucion_crout.txt\n");
            break;
        case 3:
            printf("\n--- METODO CHOLESKY ---\n");
            cholesky(N, matriz, L); // Descomposición LU
            transponer_matriz(N, L, Lt); // Transponer L
            Lx_b(N, L, b, z); // Resolución de Lz = b
            Ux_b(N, Lt, z, x); // Resolución de Ux = z
            // Comprobación Ax - b
            multiplicar_matriz_vector(N, matriz, x, Ax);
            for (int i = 0; i < N; i++) {Ax[i] -= b[i];if (fabs(Ax[i]) >= TOL) printf("La solución falló en x[%d] = %f \n", i, x[i]);}
            

            // creamos un archivo para guardar la solución
            archivo_solucion = fopen("solucion_cholesky.txt", "w");
            if (archivo_solucion == NULL) { printf("No se pudo abrir el archivo 3.\n"); return 1; }
            // escribimos la solución en el archivo
            for (int i = 0; i < N; i++) fprintf(archivo_solucion, "%f\n", x[i]);
            fclose(archivo_solucion);
            printf("\n\nLa solución fue guardada en el archivo solucion_cholesky.txt\n");
            break;
        
        default:
            printf("Saliendo del programa...\n");
            break;
        }
    } else {
        int opcion;
        printf("\nElija el método que desee aplicar:\n");    
        printf("1. Doolittle\n");
        printf("2. Crout\n");
        printf("3. Cholesky\n");
        printf("Introduzca el número de la opción: ");
        scanf("%d", &opcion);

        double **L = generar_matriz(N);
        double **U = generar_matriz(N);
        double **LU = generar_matriz(N);
        double **Lt = generar_matriz(N);
        double *x = generar_vector(N);
        double *z = generar_vector(N);
        double *Ax = generar_vector(N);
        double *comprobacion = generar_vector(N);

        switch (opcion) {
        case 1:
            printf("\n--- METODO DOOLITTLE ---\n\n");
            printf("La matriz A es:\n"); imprimir_matriz(N, matriz);
            printf("\nEl vector b es:\n"); imprimir_vector(N, b);
        
            doolittle(N, matriz, L, U); // Descomposición LU
            printf("\nLa matriz L es:\n"); imprimir_matriz(N, L);
            printf("\nLa matriz U es:\n"); imprimir_matriz(N, U);
            printf("\nLa matriz L*U es:\n"); multiplicar_matrices(N, L, U, LU); imprimir_matriz(N, LU);
            printf("\nResolucion Lz = b\n"); Lx_b(N, L, b, z); imprimir_vector(N, z); // Resolución de Lz = b
            printf("\nNuestro vector de solución es:\n"); Ux_b(N, U, z, x); imprimir_vector(N, x); // Resolución de Ux = z
            printf("\nComprobación Ax - b = 0 :\n"); 
            multiplicar_matriz_vector(N, matriz, x, Ax); 
            for (int i = 0; i < N; i++) comprobacion[i] = fabs(Ax[i] - b[i]); imprimir_vector(N, comprobacion); // Comprobación Ax - b
            for (int i = 0; i < N; i++) {if (comprobacion[i] >= TOL) printf("La solución falló en x[%d] = %f\n", i, comprobacion[i]);} // Comprobación Ax - b
            break;

        case 2:
            printf("\n--- METODO CROUT ---\n\n");
            printf("\nLa matriz A es:\n"); imprimir_matriz(N, matriz);
            printf("\nEl vector b es:\n"); imprimir_vector(N, b);

            crout(N, matriz, L, U); // Descomposición LU
            printf("\nLa matriz L es:\n"); imprimir_matriz(N, L);
            printf("\nLa matriz U es:\n"); imprimir_matriz(N, U);
            printf("\nLa matriz L*U es:\n"); multiplicar_matrices(N, L, U, LU); imprimir_matriz(N, LU);
            printf("\nResolucion Lz = b\n"); Lx_b(N, L, b, z); imprimir_vector(N, z); // Resolución de Lz = b
            printf("\nNuestro vector de solucion es:\n"); Ux_b(N, U, z, x); imprimir_vector(N, x); // Resolución de Ux = z
            printf("\nComprobacion Ax - b = 0 :\n");
            multiplicar_matriz_vector(N, matriz, x, Ax);
            for (int i = 0; i < N; i++) comprobacion[i] = fabs(Ax[i] - b[i]); imprimir_vector(N, comprobacion); // Comprobación Ax - b
            for (int i = 0; i < N; i++) {if (comprobacion[i] >= TOL) printf("La solución falló en x[%d] = %f\n", i, comprobacion[i]);} // Comprobación Ax - b
            break;
        case 3:
            printf("\n--- METODO CHOLESKY ---\n\n");
            printf("\nLa matriz A es:\n"); imprimir_matriz(N, matriz);
            printf("\nEl vector b es:\n"); imprimir_vector(N, b);

            cholesky(N, matriz, L); // Descomposición LU
            printf("\nLa matriz L es:\n"); imprimir_matriz(N, L);
            printf("\nLa matriz L^t es:\n"); transponer_matriz(N, L, Lt); imprimir_matriz(N, Lt);
            printf("\nLa resolucion Lz = b es:\n"); Lx_b(N, L, b, z); imprimir_vector(N, z); // Resolución de Lz = b
            printf("\nNuestro vector de solucion es:\n"); Ux_b(N, Lt, z, x); imprimir_vector(N, x); // Resolución de Ux = z
            printf("\nComprobacion Ax - b = 0 :\n");
            multiplicar_matriz_vector(N, matriz, x, Ax);
            for (int i = 0; i < N; i++) comprobacion[i] = fabs(Ax[i] - b[i]); imprimir_vector(N, comprobacion); // Comprobación Ax - b
            for (int i = 0; i < N; i++) {if (comprobacion[i] >= TOL) printf("La solución falló en x[%d] = %f\n", i, comprobacion[i]);} // Comprobación Ax - b
            break;
        default:
            break;
    }
    return 0;
}
}