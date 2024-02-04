#include <stdio.h>

#define INTERCAMBIAR(a, b) {double temp; temp = a; a = b; b = temp;} 
#define ABS(a) ((a) > 0 ? (a) : -(a))

int maximo_vector(int N, double v[N]);
void Ux_b(int N, double U[N][N], double b[N], double x[N]);
void gauss_pivoteo(int N, double A[N][N], double b[N], double x[N]);

int main() {
    const int N = 5;


    //                  a)  
    
     double A[5][5] = {
        {9, 7, 2, 3, 2},
        {1, 10, 9, 10, 9},
        {3, 8, 9, 1, 8},
        {3, 5, 2, 1, 4},
        {8, 8, 1, 8, 8}
    };

    double b[5] = {51, 133, 90, 43, 99};


//                       b)
/*
    double A[5][5] = {
        {2, 20, 18, 20, 18},
        {6, 16, 18, 2, 16},
        {18, 14, 4, 6, 4},
        {16, 16, 2, 16, 16},
        {6, 10, 4, 2, 8}
    };

    double b[5] = {532, 360, 204, 396, 172};
*/


    double x[N];
    printf("a)\n\n");
    printf("Matriz Original :\n");

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
    double l;
    int maximo;
    for (int k = 0; k < N - 1; k++) {
        maximo = maximo_vector(N - k, &A[k][k]) + k; 
        if (maximo != k) {
            INTERCAMBIAR(b[k], b[maximo]);
            for (int j = k; j < N; j++) INTERCAMBIAR(A[k][j], A[maximo][j]);
        }
        if (ABS(A[k][k]) < 1e-16) {
            printf("El sistema no tiene solución única\n");
            return;
        }
        for (int i = k + 1; i < N; i++) {
            l = A[i][k] / A[k][k];
            for (int j = k + 1; j < N; j++) A[i][j] = A[i][j] - l * A[k][j];
            b[i] = b[i] - l * b[k];
        }
    }
    Ux_b(N, A, b, x);
}
