#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

void gradiente_conjugado(int N, double tolerancia, double **matriz, double *b, double *vector);

int main(){
    int N = 3;
    double **matriz = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) matriz[i] = (double *)malloc(N * sizeof(double));

    matriz[0][0] = 3.0; matriz[0][1] = 1.0; matriz[0][2] = 4.0;
    matriz[1][0] = 1.0; matriz[1][1] = 5.0; matriz[1][2] = 9.0;
    matriz[2][0] = 4.0; matriz[2][1] = 9.0; matriz[2][2] = 7.0;

    double *b = (double *)malloc(N * sizeof(double));


    b[0] = 3.0; b[1] = 1.0; b[2] = 2.0;

    double *vector = (double *)malloc(N * sizeof(double));
    vector[0] = 1.0; vector[1] = 1.0; vector[2] = 1.0;
    // vector[0] = 2.0; vector[1] = 1.0;

    // imprimimos el sistema de ecuaciones
    printf("Sistema de ecuaciones:\n");
    imprimir_matriz(N, matriz);
    printf("\n");
    printf("Vector b:\n");
    imprimir_vector(N, b);
    printf("\n");
    printf("Vector inicial:\n");
    imprimir_vector(N, vector);
    printf("\n");

    // hacemos el método de gradiente conjugado
    double tolerancia = 1e-16;
    gradiente_conjugado(N, tolerancia, matriz, b, vector);

    // imprimimos la solución
    printf("\nSolución:\n");
    imprimir_vector(N, vector);

    // liberamos memoria 
    liberar_matriz(N, matriz);
    free(b);
    free(vector);

    return 0;
}

void gradiente_conjugado(int N, double tolerancia, double **matriz, double *b, double *vector) {
    int iteraciones = 2 * N;
    double *Ax = (double *)malloc(N * sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *r_nuevo = (double *)malloc(N * sizeof(double));
    double *p = (double *)malloc(N * sizeof(double));
    double *paso_Apk = (double *)malloc(N * sizeof(double));

    // Inicializa rk, pk y rk_nuevo a cero
    double ak, bk;

    // Calcula Ax
    multiplicar_matriz_vector(N, matriz, vector, Ax);

    // Calcula r0
    for (int i = 0; i < N; i++) r[i] = b[i] - Ax[i];

    // Calcula la norma de r0
    double max = 0.0;
    for (int i = 0; i < N; i++) if (fabs(r[i]) > max) max = fabs(r[i]);

    // Comprueba la condición de convergencia
    if (max < tolerancia) {
        printf("Gradiente Conjugado: El vector inicial es una solución aproximada\n");
        free(Ax);
        free(r);
        free(r_nuevo);
        free(p);
        free(paso_Apk);
        return;
    }

    // Inicializa p0
    for (int i = 0; i < N; i++) p[i] = r[i];

    for (int iter = 0; iter < iteraciones; iter++) {
        // Calcula ak
        multiplicar_matriz_vector(N, matriz, p, paso_Apk);
        ak = producto_punto(N, r, r) / producto_punto(N, p, paso_Apk);

        // Calcula x_k+1
        for (int i = 0; i < N; i++) vector[i] += ak * p[i];

        // Calcula r_k+1
        for (int i = 0; i < N; i++) r_nuevo[i] = r[i] - ak * paso_Apk[i];

        // Comprueba la condición de convergencia
        if (sqrt(producto_punto(N, r_nuevo, r_nuevo)) < tolerancia) {
            printf("Gradiente conjugado: el método convergió después de %d iteraciones.\n", iter + 1);
            free(Ax);
            free(r);
            free(r_nuevo);
            free(p);
            free(paso_Apk);
            return;
        }

        // Calcula bk
        bk = producto_punto(N, r_nuevo, r_nuevo) / producto_punto(N, r, r);

        // Actualiza los valores para la siguiente iteración
        for (int i = 0; i < N; i++) p[i] = r_nuevo[i] + bk * p[i];

        // Actualiza r
        for (int i = 0; i < N; i++) r[i] = r_nuevo[i];

    }

    printf("Gradiente conjugado: El método no convergió después de %d iteraciones.\n", 2 * iteraciones);

    // Liberar memoria
    free(Ax);
    free(r);
    free(r_nuevo);
    free(p);
    free(paso_Apk);
}