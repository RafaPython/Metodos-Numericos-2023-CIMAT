#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

/*

    * Método de factorización QR

*/

// void QR(int N, int M, double **matriz, double **Q, double **R);
void comprobacion_QR(int N,int M, double **matriz, double **Q, double **R, double tolerancia);

int main(){

    int N = 3;
    int M = 3;

    // ejemplo de una matriz 
    double **matriz = (double **)malloc(N *sizeof(double *));
    for(int i = 0; i < N; i++) matriz[i] = (double*)malloc(N*sizeof(double));

    matriz[0][0] = 12.0; matriz[0][1] = -51.0; matriz[0][2] = 4.0;
    matriz[1][0] = 6.0; matriz[1][1] = 167.0; matriz[1][2] = -68.0;
    matriz[2][0] = -4.0; matriz[2][1] = 24.0; matriz[2][2] = -41.0;

    // abrimos memoria para las matrices Q y R 

    double **Q = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) Q[i] = (double *)malloc(N * sizeof(double));

    double **R = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) R[i] = (double *)malloc(M * sizeof(double));

    // realizamos la factorización QR
    QR(N, M, matriz, Q, R);

    // imprimimos las matrices Q y R
    printf("Matriz Q:\n");
    imprimir_matriz(N, Q);
    printf("\n");

    // limpiamos el buffer
    
    // printf("Matriz R:\n");
    // imprimir_matriz(N, R);
    // printf("\n");
    //imprimimo R
    
    // imprimimos la matriz original 
    printf("Matriz original:\n");
    imprimir_matriz(N, matriz);
    printf("\n");

    // comprobamos que la factorización es correcta

    comprobacion_QR(3, 3, matriz, Q, R, 1e-10);

    // liberamos memoria
    liberar_matriz(N, matriz);
    liberar_matriz(N, Q);
    liberar_matriz(N, R);
    


    return 0;
}

// void QR(int N, int M, double **matriz, double **Q, double **R) {

//     // Copiar la matriz original en una nueva matriz
//     double **A = (double **)malloc(N * sizeof(double *));
//     for(int i = 0; i < N; i++) A[i] = (double *)malloc(M * sizeof(double));
//     for(int i = 0; i < N; i++) for(int j = 0; j < M; j++) A[i][j] = matriz[i][j];

//     double norma;
//     for (int i = 0; i < N; i++) {
//         // Calcular la norma de la columna i
//         norma = 0.0;
//         for (int j = 0; j < M; j++) norma += A[j][i] * A[j][i];
//         norma = sqrt(norma);

//         // Actualizar la columna i de R y la columna i de Q
//         R[i][i] = norma;
//         for (int j = 0; j < M; j++) Q[j][i] = A[j][i] / R[i][i];

//         // Actualizar el resto de R y Q
//         for (int j = i + 1; j < N; j++) {
//             R[i][j] = 0.0;
//             for (int k = 0; k < M; k++) R[i][j] += Q[k][i] * A[k][j];
//             for (int k = 0; k < M; k++) A[k][j] -= R[i][j] * Q[k][i];
//         }
//     }

//     liberar_matriz(N, A);
// }

void comprobacion_QR(int N,int M, double **matriz, double **Q, double **R, double tolerancia){
    double **aux1 = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) aux1[i] = (double *)malloc(N* sizeof(double));




    multiplicar_matrices(N, Q, R, aux1);
    // imprimir_matriz(N, aux1);

    int contador = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(fabs(aux1[i][j] - matriz[i][j]) > tolerancia){
                printf("La factorización QR no es correcta en  A[%d][%d]\n",i,j);
                contador++;
            }
        }
    }

    if(contador == 0) printf("La factorización QR es correcta.\n");
    liberar_matriz(N, aux1);

}
