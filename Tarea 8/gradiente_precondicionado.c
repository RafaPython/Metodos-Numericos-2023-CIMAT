#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

void gradiente_conjugado_precondicionado(int N, double tolerancia, double **matriz, double *b, double *vector, double **precondicionador);

int main(){
    int N = 3;
    double **matriz = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) matriz[i] = (double *)malloc(N * sizeof(double));

    matriz[0][0] = 3.0; matriz[0][1] = 1.0; matriz[0][2] = 6.0;
    matriz[1][0] = 1.0; matriz[1][1] = 5.0; matriz[1][2] = 9.0;
    matriz[2][0] = 6.0; matriz[2][1] = 9.0; matriz[2][2] = 7.0;

    double *b = (double *)malloc(N * sizeof(double));
    b[0] = 3.0; b[1] = 1.0; b[2] = 2.0;

    double *vector = (double *)malloc(N * sizeof(double));
    vector[0] = 1.0; vector[1] = 1.0; vector[2] = 1.0;

    double **precondicionador = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) precondicionador[i] = (double *)malloc(N * sizeof(double));

        // calculamos el precondicionador M ^-1 = D^-1
    
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) precondicionador[i][j] = (i == j) ? 1.0 / matriz[i][j] : 0.0;


    // Inicializa la matriz de precondicionamiento (puede ser una matriz diagonal)
    // for(int i = 0; i < N; i++) {
    //     for(int j = 0; j < N; j++) {
    //         if(i == j)
    //             precondicionador[i][j] = 1.0;  // Usando una matriz identidad como precondicionador por defecto
    //         else
    //             precondicionador[i][j] = 0.0;
    //     }
    // }

    // Imprimimos el sistema de ecuaciones
    printf("Sistema de ecuaciones:\n");
    imprimir_matriz(N, matriz);
    printf("\n");
    printf("Vector b:\n");
    imprimir_vector(N, b);
    printf("\n");
    printf("Vector inicial:\n");
    imprimir_vector(N, vector);
    printf("\n");

    // Hacemos el método de gradiente conjugado precondicionado
    double tolerancia = 1e-16;


    gradiente_conjugado_precondicionado(N, tolerancia, matriz, b, vector, precondicionador);

    // Imprimimos la solución
    printf("\nSolución:\n");
    imprimir_vector(N, vector);

    // Liberamos memoria
    liberar_matriz(N, matriz);
    free(b);
    free(vector);
    for(int i = 0; i < N; i++) free(precondicionador[i]);
    free(precondicionador);

    return 0;
}

void gradiente_conjugado_precondicionado(int N, double tolerancia, double **matriz, double *b, double *vector, double **precondicionador) {
    int iteraciones = 2 * N;
    double *Ax = (double *)malloc(N * sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *r_nuevo = (double *)malloc(N * sizeof(double));
    double *z = (double *)malloc(N * sizeof(double));
    double *p = (double *)malloc(N * sizeof(double));
    double *paso_Apk = (double *)malloc(N * sizeof(double));

    // Inicializa rk, pk y rk_nuevo a cero
    double ak, bk;

    // Calcula Az
    multiplicar_matriz_vector(N, matriz, vector, Ax);

    // Calcula r0
    for (int i = 0; i < N; i++) r[i] = b[i] - Ax[i];

    // Aplica precondicionamiento: z = M^(-1) * r
    for (int i = 0; i < N; i++) {
        z[i] = r[i] / precondicionador[i][i];
    }

    // Calcula la norma de r0
    double max = 0.0;
    for (int i = 0; i < N; i++) if (fabs(r[i]) > max) max = fabs(r[i]);

    // Comprueba la condición de convergencia
    if (max < tolerancia) {
        printf("Gradiente Conjugado Precondicionado: El vector inicial es una solución aproximada\n");
        free(Ax);
        free(r);
        free(r_nuevo);
        free(z);
        free(p);
        free(paso_Apk);
        return;
    }

    // Inicializa p0
    for (int i = 0; i < N; i++) p[i] = z[i];

    for (int iter = 0; iter < iteraciones; iter++) {
        // Calcula ak
        multiplicar_matriz_vector(N, matriz, p, paso_Apk);

        double rz_inner = producto_punto(N, z, r);
        ak = rz_inner / producto_punto(N, p, paso_Apk);

        // Calcula x_k+1
        for (int i = 0; i < N; i++) vector[i] += ak * p[i];

        // Calcula r_k+1
        for (int i = 0; i < N; i++) r_nuevo[i] = r[i] - ak * paso_Apk[i];

        // Comprueba la condición de convergencia
        if (sqrt(producto_punto(N, r_nuevo, r_nuevo)) < tolerancia) {
            printf("Gradiente conjugado precondicionado: el método convergió después de %d iteraciones.\n", iter + 1);
            free(Ax);
            free(r);
            free(r_nuevo);
            free(z);
            free(p);
            free(paso_Apk);
            return;
        }

        // Aplica precondicionamiento: z_k+1 = M^(-1) * r_k+1
        for (int i = 0; i < N; i++) {
            z[i] = r_nuevo[i] / precondicionador[i][i];
        }

        // Calcula bk
        double rz_nuevo_inner = producto_punto(N, z, r_nuevo);
        bk = rz_nuevo_inner / rz_inner;

        // Actualiza los valores para la siguiente iteración
        for (int i = 0; i < N; i++) p[i] = z[i] + bk * p[i];

        // Actualiza r
        for (int i = 0; i < N; i++) r[i] = r_nuevo[i];
    }

    printf("Gradiente conjugado precondicionado: El método no convergió después de %d iteraciones.\n", 2 * iteraciones);

    // Liberar memoria
    free(Ax);
    free(r);
    free(r_nuevo);
    free(z);
    free(p);
    free(paso_Apk);
}
