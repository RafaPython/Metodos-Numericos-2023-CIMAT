#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

// void rayleigh(int N, double sigma, double **matriz, double *v0, double tolerancia, int iteraciones);

int main(){
    int N = 3;
    double **matriz = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) matriz[i] = (double *)malloc(N * sizeof(double));

    matriz[0][0] = 3.0; matriz[0][1] = -0.1; matriz[0][2] = -0.2;
    matriz[1][0] = -0.1; matriz[1][1] = 7.0; matriz[1][2] =-0.3;
    matriz[2][0] = -0.2; matriz[2][1] = -0.3; matriz[2][2] = 10.0;

    // inicizamos el vector inicial
    double *v0 = (double *)malloc(N * sizeof(double));
    v0[0] = 1.0; v0[1] = 1.0; v0[2] = 1.0;

    double sigma = 1.5 ; 
    double tolerancia = 1e-16;
    int iteraciones = 100;

    // imprimimos la matriz original
    printf("Matriz original:\n");
    imprimir_matriz(N, matriz);
    printf("\n");

    printf("Vector inicial:\n");
    imprimir_vector(N, v0);
    printf("\n");

    printf("Sigma original: %lf\n", sigma);
    printf("\n");
    // hacemos el método de Rayleig h
    
    sigma = rayleigh(N, sigma, matriz, v0, tolerancia, iteraciones);

    printf("Vector final:\n");
    imprimir_vector(N, v0);

    printf("Sigma final: %lf\n", sigma);

    // liberamos memoria
    liberar_matriz(N, matriz);
    free(v0);



    return 0;
}


// void rayleigh(int N, double sigma, double **matriz, double *v0, double tolerancia, int iteraciones){

//     // abirmos espacio para la matriz B 
//     double **B = (double **)malloc(N * sizeof(double *));
//     for(int i = 0; i < N; i++) B[i] = (double *)malloc(N * sizeof(double));

//     // abrimos espacio para las matrices B = LU 
//     double **U = (double **)malloc(N * sizeof(double *));
//     for(int i = 0; i < N; i++) U[i] = (double *)malloc(N * sizeof(double));
    
//     double **L = (double **)malloc(N * sizeof(double *));
//     for(int i = 0; i < N; i++) L[i] = (double *)malloc(N * sizeof(double));

//     double *Ax = (double *)malloc(N * sizeof(double *));

//     double *v1 = (double *)malloc(N * sizeof(double));
//     for(int i = 0; i < N; i++) v1[i] = 1.0;
//     double *y = (double *)malloc(N * sizeof(double));
//     double *aux = (double *)malloc(N * sizeof(double));

//     double norma;

//     for (int i = 1; i <= iteraciones; i++){

//         // realizamos la resta de la matriz A - sigma * I
//         for(int i = 0; i < N; i++)
//             for (int j = 0; j < N; j++) 
//                 B[i][j] = (i == j) ?  matriz[i][j] - sigma : matriz[i][j]; 

//         // ahora resolvemos el sistema de ecuaciones 
//         doolittle(N, B, L, U);

//         Lx_b(N, L, v0, y);
//         Ux_b(N, U, y, v1);        

//         norma = sqrt(producto_punto(N, v1, v1));

//         for(int i = 0; i < N; i++) v1[i] = v1[i] / norma;

//         // calculamos abs(v1 - v0)
        
//         for(int i = 0; i < N; i++) aux[i] = v1[i] - v0[i];

//         norma = sqrt(producto_punto(N, aux, aux));

//         if(norma < tolerancia){
//             printf("El método de Rayleigh convergió en %d iteraciones.\n", i);
//             printf("El autovalor es: %lf\n", sigma);
//             printf("El autovector es:\n");
//             imprimir_vector(N, v1);
//             printf("\n");

//             // liberamos la memoria
//             liberar_matriz(N, B);
//             liberar_matriz(N, L);
//             liberar_matriz(N, U);
//             free(Ax);
//             free(v1);
//             free(y);
//             free(aux);

//             return;
//         }

//         // actualizamos el vector v0 = v1
//         for(int i = 0; i < N; i++) v0[i] = v1[i];

//         // actualizamos el sigma
//         multiplicar_matriz_vector(N, matriz, v0, Ax);
    
//         sigma = producto_punto(N, v0, Ax) ;
//     }

//     printf("El método de Rayleigh no convergió en %d iteraciones.\n", iteraciones);

//     // liberar memoria
//     liberar_matriz(N, B);
//     liberar_matriz(N, L);
//     liberar_matriz(N, U);
//     free(Ax);
//     free(v1);
//     free(y);
//     free(aux);
// }