/*
Problema 3 Tarea 3 Métodos Numéricos 

Autor: Rafael Alejandro García Ramírez
email: rafael.ramirez@cimat.mx

*/

#include <stdio.h>

void construir_matriz_2triangular_superior(int N, double A[N][N], double b[N]);
void Ux_b(int N, double U[N][N], double b[N] ,double x[N]);

int main() {
    // Ejemplo
    const int N = 3;

    double A[3][3] = {{3, 1, 3}, {2, 4, 3}, {7, 3, 5}};
    double b[3] = {2, 5, 8};

    double x[N];
    printf("Matriz Original:\n");
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%f ", A[i][j]);
        printf("\n");
    }

    construir_matriz_2triangular_superior(N, A, b);

    printf("Matriz triangular superior:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%f ", A[i][j]);
        printf("\n");
    }

    Ux_b(N, A, b, x);
    printf("Soluciones:\n");
    for (int i = 0; i < N; i++) printf("x[%d] = %f\n", i, x[i]);
    
    return 0;
}

void construir_matriz_2triangular_superior(int N, double A[N][N], double b[N]) {
    /* Construye una matriz A y un vector b en una forma triangular superior.
     *
     * @ N       El tamano de la matriz cuadrada A y el vector b.
     * @ A       La matriz cuadrada de tamano N x N.
     * @ b       El vector de terminos independientes de tamano N.
     *
     * Esta funcion modifica la matriz A y el vector b para llevarlos a una forma triangular superior.
     * Puede utilizarse como paso intermedio para resolver sistemas de ecuaciones lineales.
     */
    
    double l;
    for (int k = 0; k < N; k++) {
        for (int i = k + 1; i < N; i++) {
            l = A[i][k] / A[k][k];
            for (int j = k; j < N; j++) A[i][j] = A[i][j] - l * A[k][j];
            b[i] = b[i] - l * b[k];
        }
    }
}

void Ux_b(int N, double U[N][N], double b[N] ,double x[N]){
    /* Resuelve un sistema de ecuaciones lineales triangular superior Ux = b para x.
     *
     * @ N       El tamano de la matriz triangular superior U y los vectores b y x.
     * @ U       La matriz triangular superior de tamano N x N.
     * @ b       El vector de terminos independientes de tamano N.
     * @ x       El vector de solucion, que se actualizara con la solucion despues de la llamada.
     *
     * Esta funcion utiliza sustitucion hacia atras para resolver el sistema de ecuaciones lineales.
     * La matriz U debe ser triangular superior, y el vector x debe contener una estimacion inicial
     * de la solucion, que se actualizara con la solucion final.
     */

    double acumulador;
    for(int i = N - 1; i >= 0; i--){
        acumulador = 0;
        for (int k =  N; k > i; k--) acumulador += U[i][k]*x[k];
        x[i] = (b[i] - acumulador) / U[i][i];
    }  
}