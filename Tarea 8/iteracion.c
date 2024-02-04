#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"
// void restar_proyeccion(int p, int q, double **eigenvectores, double*vector);
// void iteracion_subespacio(int N, double **matriz,double **eigenvalores, double **eigenvectores, double tolerancia, int iteraciones);

int main(){

    int N = 3;
    double **matriz = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) matriz[i] = (double *)malloc(N * sizeof(double));

    matriz[0][0] = 3.0; matriz[0][1] = -0.1; matriz[0][2] = -0.2;
    matriz[1][0] = -0.1; matriz[1][1] = 7.0; matriz[1][2] =-0.3;
    matriz[2][0] = -0.2; matriz[2][1] = -0.3; matriz[2][2] = 10.0;

    double tolerancia = 1e-10;
    int iteraciones = 10000;
    
    printf("Matriz original:\n");
    imprimir_matriz(N, matriz);
    printf("\n");

    double **eigenvalores = generar_matriz(N);  
    double **eigenvectores = generar_matriz(N);
    // inicializamos los eigenvectores como la matriz identidad
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) eigenvectores[i][j] = (i == j) ? 1.0 : 0.0;

    iteracion_subespacio(N, matriz, eigenvalores, eigenvectores, tolerancia, iteraciones);

    printf("Eigenvalores:\n");
    imprimir_matriz(N, eigenvalores);
    printf("\n");

    printf("Eigenvectores:\n");
    imprimir_matriz(N, eigenvectores);
    printf("\n");


    // liberamos memoria
    liberar_matriz(N, matriz);
    liberar_matriz(N, eigenvalores);
    liberar_matriz(N, eigenvectores);


    return 0;
}

// void restar_proyeccion(int p, int q, double **eigenvectores, double*vector){
//     double dummy;
//     for(int i = 0; i < q; i++){
//         dummy = 0;
//         for(int j = 0; j < p; j++) dummy += eigenvectores[j][i]*vector[j];
//         for(int j = 0; j < p; j++) vector[j] -= dummy*eigenvectores[j][i];
//     }
// }
// void normalizar(double *vector, int longitud) {
//     double sumaCuadrados = 0;
    
//     // Calcula la suma de los cuadrados de los elementos del vector
//     for (int i = 0; i < longitud; i++) {
//         sumaCuadrados += vector[i] * vector[i];
//     }
    
//     // Calcula la magnitud
//     double magnitud = sqrt(sumaCuadrados);
    
//     // Normaliza el vector dividiendo cada elemento por la magnitud
//     for (int i = 0; i < longitud; i++) {
//         vector[i] /= magnitud;
//     }
// }
// void iteracion_subespacio(int N, double **matriz,double **eigenvalores, double **eigenvectores, double tolerancia, int iteraciones){
//     double **p = generar_matriz(N);
//     double **pt = generar_matriz(N);
    
//     double **L = generar_matriz(N);
//     double **U = generar_matriz(N);
//     double **aux = generar_matriz(N);

//     double *v0 = (double *)malloc(N*sizeof(double));
//     double *v1 = (double *)malloc(N*sizeof(double));
//     double *y = (double *)malloc(N*sizeof(double));

//     // factorizamos A = LU
//     doolittle(N, matriz, L, U);

//     int diagonal, identidad;

//     // inicializamos la matriz p
//     for(int i = 0; i < N; i++) p[i][i] = 1.0;

//     for(int iter = 0; iter < iteraciones; iter++){

//         // Realizamos potencia inversa
//         for(int i = 0; i < N; i++){ 
//             for(int j = 0; j < N; j++) v0[j] = p[j][i];
//             Lx_b(N, L, v0, y);
//             Ux_b(N, U, y, v1);
//             restar_proyeccion(N, i, p, v1);
//             normalizar(v1, N);
//             for(int j = 0; j < N; j++) p[j][i] = v1[j]; 

//         }

//         trasponer_matriz(N,p,pt);
        
//         multiplicar_matrices(N, matriz, p, aux);
//         multiplicar_matrices(N,pt , aux, eigenvalores);

//         // verificamos si eigenvalores es diagonal 
//         diagonal = 1;
//         for (int j = 0; j < N; j++) {
//             for (int k = 0; k < N; k++) {
//                 if (j != k && fabs(eigenvalores[j][k]) > tolerancia) {
//                     diagonal = 0;  // eigenvalores no es diagonal
//                 }
//             }
//         }
//         if (diagonal == 1) break;

//         jacobi_eigen(N,eigenvalores, eigenvectores, tolerancia,iteraciones);

//         identidad = 1;
//         for (int j = 0; j < N; j++) {
//             for (int k = 0; k < N; k++) {
//                 if (j == k) {
//                     // Verificar si el elemento en la diagonal principal es 1
//                     if (fabs(eigenvectores[j][k] - 1.0) > tolerancia) {
//                         identidad = 0;  // eigenvectores no es la matriz identidad
//                     }
//                 } else {
//                     // Verificar si los elementos fuera de la diagonal son 0
//                     if (fabs(eigenvectores[j][k]) > tolerancia) {
//                         identidad =  0;  // eigenvectores no es la matriz identidad
//                     }
//                 }
//             }
//         }

//         if (identidad == 1) break;

//         if(iter == iteraciones - 1){
//             printf("Error Iteracion subespacio: el método falló después de %d iteraciones",iteraciones);
//             liberar_matriz(N, p);
//             liberar_matriz(N, pt);
//             liberar_matriz(N, L);
//             liberar_matriz(N, U);
//             liberar_matriz(N, aux);
//             free(v0);
//             free(v1);
//             free(y);

//             return;
//         }

//         multiplicar_matrices(N, matriz, p, aux);


//         for(int j = 0; j < N; j++) for(int k = 0; k < N; k++) p[j][k] = aux[j][k];
        
//     }

//     for(int i = 0; i < N; i++){
//         for(int j = 0; j < N; j++){
//             if (i != j){
//                 if (fabs(eigenvalores[i][j]) <= tolerancia) eigenvalores[i][j] = 0.0;
//             }
//         }
//     }

//     liberar_matriz(N, p);
//     liberar_matriz(N, pt);
//     liberar_matriz(N, L);
//     liberar_matriz(N, U);
//     liberar_matriz(N, aux);
//     free(v0);
//     free(v1);
//     free(y);
// }


