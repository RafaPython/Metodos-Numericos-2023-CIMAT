#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "solvers.h"

void matriz_pentadiagonal(int N, double **matriz, double *b);

void imprimir_matriz_n(int N, int M, double **matriz){
    // funcion que imprime una sección de nxn de una matriz más grande
    for(int i = N; i < M; i++){
        for (int j = N; j < M; j++){
            printf("%.2f ", matriz[i][j]);
        }
        printf("\n");
    }
}
void imprimir_vector_n(int N, int M, double *vector){
    // funcion que imprime una sección de nxn de una matriz más grande
    for(int i = N; i < M; i++) printf("%.6f ", vector[i]);
}

int diferencia_vectores(int N, double tolerancia, double *vector1, double *vector2){
    // funcion que devuelve la diferencia entre dos vectores
    int contador = 0;
    for(int i = 0; i < N; i++) if(fabs(vector1[i] - vector2[i]) > tolerancia) contador++;
    return contador;
}

int main(){


    int N = 2000;
    double **matriz = generar_matriz(N);
    double *b = generar_vector(N);

    matriz_pentadiagonal(N, matriz, b);

    double tol = 1e-4;
    int iteracionesMAX = 100000;

    printf("Tolerancia --> %.16f\n", tol);
    printf("Iteraciones MAX --> %d\n", iteracionesMAX);

    
    // void gauss_seidel(int N, double TOL, double **matriz, double *b, double *x0, int iteracionesMAX);

    double *x0 = generar_vector(N);
    for(int i = 0; i < N; i++) x0[i] = 1.0;

//     ---- GAUSS SEIDEL ----
///////////////////////////////////////////////////
    printf("\n\nGAUSS SEIDEL\n\n");
    gauss_seidel(N, tol, matriz, b, x0, iteracionesMAX);
    

    // imprimimos algunos valores de la solución
    printf("\nPrimeras 23 iteraciones:\n");
    imprimir_vector_n(0, 23, x0);


    // imprimimos algunos valores de la solución
    printf("\n\nUltimas 23 iteraciones:\n");
    imprimir_vector_n(N - 23, N, x0);

    // comprobamos que la solución es correcta
    printf("\n\nLa cantidad de elementos diferentes de cero en Ax - b = 0 es de: %d\n", comprobar_solucion(N, matriz, x0, b, tol));


//////////////////////////////////////////////////////////////////////

    printf("\n\nGRADIENTE CONJUGADO PRECONDICIONADO\n\n");
    double *x1 = generar_vector(N);
    for(int i = 0; i < N; i++) x1[i] = rand() % 2;

    double **precondcionador = generar_matriz(N);
    for (int i = 0; i < N; i++) for(int j = 0; j < N; j++) precondcionador[i][j] = (i == j) ? matriz[i][j] : 0.0;

    // void gradiente_conjugado_precondicionado(int N, double tolerancia, double **matriz, double *b, double *vector, double **precondicionador)
    gradiente_conjugado_precondicionado(N, tol, matriz, b, x1, precondcionador);

    // imprimimos algunos valores de la solución
    printf("\nPrimeras 23 iteraciones:\n");
    imprimir_vector_n(0, 23, x1);


    // imprimimos algunos valores de la solución
    printf("\n\nUltimas 23 iteraciones:\n");
    imprimir_vector_n(N - 23, N, x1);


    // comprobamos que la solución es correcta
    printf("\n\nLa cantidad de elementos diferentes de cero en Ax - b = 0 es de: %d\n", comprobar_solucion(N, matriz, x0, b, tol));

    // comparamos ambos resultados
    printf("\n\n\nLa cantidad de elementos diferentes de cero en x0 - x1 = 0 es de: %d\n", diferencia_vectores(N, tol, x0, x1));


    ///////////////////////////////////////////////


    // potencia_k_valores_mayores(int N, int k, double **matriz, double **vectores, double *eigenvalores, double tolerancia, int iteraciones)

    // printf("\n\nPOTENCIA INVERSA\n\n");
    // int k = 10;
    // double **vectores = generar_matriz(N);
    // double *eigenvalores = generar_vector(N);
    
    // potencia_k_valores_mayores(N, k, matriz, vectores, eigenvalores, tol, iteracionesMAX);

    // // imprimimos algunos valores de la solución
    // printf("\nPrimeros 10 eigenvalores:\n");
    // imprimir_vector_n(0, 10, eigenvalores);

    // printf("\n\nEigenvectores: \n");
    // imprimir_matriz_n(0, 10, vectores);


    // void jacobi_eigen(int N, double **matriz, double **eigenvectores, double tolerancia, int max_iteraciones)

    // printf("\n\nJACOBI EIGEN\n\n");

    // double **eigenvectores = generar_matriz(N);
    // jacobi_eigen(N, matriz, eigenvectores, tol, iteracionesMAX);

    // printf("\n\nPrimeros 10 eigenvectores: \n");
    // imprimir_matriz_n(0, 10, eigenvectores);

    // printf("\n\nPrimeros 10 eigenvalores:\n");
    // for(int i = 0; i < 10; i++) printf("%.16f ", matriz[i][i]);

    // void iteracion_subespacio(int N, double **matriz,double **eigenvalores, double **eigenvectores, double tolerancia, int iteraciones)

    // double **eigenvectores = generar_matriz(N);
    // double **eigenvalores = generar_matriz(N);
    // printf("\n\nITERACION SUBESPACIO\n\n");
    // iteracion_subespacio(N, matriz, eigenvalores, eigenvectores, tol, iteracionesMAX);

    // printf("\n\nPrimeros 10 eigenvectores: \n");
    // imprimir_matriz_n(0, 10, eigenvectores);
    
    // printf("\n\nPrimeros 10 eigenvalores:\n");
    // for(int i = 0; i < 10; i++) printf("%.16f ", matriz[i][i]);

    ///////////////////////////////////////////
    


    liberar_matriz(N, matriz);
    liberar_matriz(N, precondcionador);
    // liberar_matriz(N, vectores);
    // free(eigenvalores);
    // liberar_matriz(N, eigenvectores);
    // liberar_matriz(N, eigenvalores);
    free(b);
    free(x0);
    free(x1);





}

void matriz_pentadiagonal(int N, double **matriz, double *b){
    /*
    Función que genera una matriz pentadiagonal de 2000x2000 :
    con el esquema −4x_i−2 − 8x_i−1 + 40x_i − 8x_i+1 − 4x_i+2 = 100, para i = 3, ..., 1998 
    Los demás elementos para i = 1, 2, 1999, 2000 son 0
    */
   // rellenamos de ceros la matriz 
   for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) matriz[i][j] = 0;
   
   for (int i = 2; i < N - 2; i++){
    matriz[i][i - 2] = -4;
    matriz[i][i - 1] = -8;
    matriz[i][i] = 40;
    matriz[i][i + 1] = -8;
    matriz[i][i + 2] = -4;
   }
    
    // generamos el vector b
    for (int i = 2; i < N - 1; i++) b[i] = 100;

    // sustituimos la primera ecuación 40x1 − 8x2 − 4x3 = 20
    matriz[0][0] = 40;
    matriz[0][1] = -8;
    matriz[0][2] = -4;
    b[0] = 20;

    // segunda ecuacion −8x1 + 40x2 − 8x3 − 4x4 = 50
    matriz[1][0] = -8;
    matriz[1][1] = 40;
    matriz[1][2] = -8;
    matriz[1][3] = -4;
    b[1] = 50;

    // la penúltima ecuación −4x1997 − 8x1998 + 40x1999 − 8x2000 = 50
    matriz[N - 2][N - 4] = -4;
    matriz[N - 2][N - 3] = -8;
    matriz[N - 2][N - 2] = 40;
    matriz[N - 2][N - 1] = -8;
    b[N - 2] = 50;

    // ultima ecuación −4x1998 − 8x1999 + 40x2000 = 20
    matriz[N - 1][N - 3] = -4;
    matriz[N - 1][N - 2] = -8;
    matriz[N - 1][N - 1] = 40;
    b[N - 1] = 20;
    
}