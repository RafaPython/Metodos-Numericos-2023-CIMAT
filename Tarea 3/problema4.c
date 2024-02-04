/*
Problema 4 Tarea 3 Métodos Numéricos 

Autor: Rafael Alejandro García Ramírez
email: rafael.ramirez@cimat.mx

*/

#include <stdio.h>

#define INTERCAMBIAR(a, b) {double temp; temp = a; a = b; b = temp;} 
#define ABS(a) ((a) > 0 ? (a) : -(a))

int maximo_vector(int N, double v[N]);
void Ux_b(int N, double U[N][N], double b[N], double x[N]);
void gauss_pivoteo(int N, double A[N][N], double b[N], double x[N]);

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
    printf("Vector b:\n");
    for (int i = 0; i < N; i++) printf("%f\n", b[i]);

    gauss_pivoteo(N, A, b, x);
    printf("Soluciones:\n");
    for (int i = 0; i < N; i++) printf("x[%d] = %f\n", i, x[i]);
    return 0;
}

int maximo_vector(int N, double v[N]) {
    // Encuentra el maximo de un vector de tamano N. Retorna el indice del maximo.
    int max = 0;
    for (int i = 1; i < N; i++) max = (v[i] > v[max]) ? i : max;
    return max;
}

void Ux_b(int N, double U[N][N], double b[N], double x[N]) {
    double acumulador;
    for (int i = N - 1; i >= 0; i--) {
        acumulador = 0;
        for (int k = N - 1; k > i; k--) acumulador += U[i][k] * x[k];
        x[i] = (b[i] - acumulador) / U[i][i];
    }
}

void gauss_pivoteo(int N, double A[N][N], double b[N], double x[N]) {
    /* Resuelve un sistema de ecuaciones lineales utilizando el método de eliminación de Gauss con pivoteo parcial.
     *
     * @ N       El tamano de la matriz cuadrada A, el vector b y el vector de solucion x.
     * @ A       La matriz cuadrada de coeficientes de tamano N x N.
     * @ b       El vector de terminos independientes de tamano N.
     * @ x       El vector de solucion, que se actualizara con la solucion despues de la llamada.
     *
     * Esta funcion resuelve el sistema de ecuaciones lineales Ax = b utilizando eliminacion de Gauss
     * con pivoteo parcial para mejorar la estabilidad numerica.
     */

    double l; int maximo;
    for (int k = 0; k < N - 1; k++) {
        // Intercambia filas para llevar el elemento de mayor magnitud a la diagonal
        maximo = maximo_vector(N - k, &A[k][k]) + k; 
        if (maximo != k) {
            INTERCAMBIAR(b[k], b[maximo]);
            for (int j = k; j < N; j++) INTERCAMBIAR(A[k][j], A[maximo][j]);
        }
        // Si el elemento de la diagonal es cero, el sistema no tiene solucion unica
        if (ABS(A[k][k]) < 1e-16) {
            printf("El sistema no tiene solucion unica\n");
            return;
        }
        // Eliminacion de Gauss
        for (int i = k + 1; i < N; i++) {
            l = A[i][k] / A[k][k];
            for (int j = k + 1; j < N; j++) A[i][j] = A[i][j] - l * A[k][j];
            b[i] = b[i] - l * b[k];
        }
    }
    Ux_b(N, A, b, x);
}