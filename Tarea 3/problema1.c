/*
Problema 1 Tarea 3 Métodos Numéricos 

Autor: Rafael Alejandro García Ramírez
email: rafael.ramirez@cimat.mx

*/


#include <stdio.h>

void Lx_b(int N, double L[N][N], double b[N] ,double x[N]);

int main(){
    // Ejemplo
    const int N = 3;

    double L[3][3] = {{3,0,0}, {2,4,0}, {1,3,5}};
    double b[3] = {2,5,8};

    double x[N];
    Lx_b(N, L, b, x);
    for (int i = 0; i < N; i++) printf("x[%d] = %f\n", i, x[i]);
    return 0;
}

void Lx_b(int N, double L[N][N], double b[N], double x[N]) {
    /* Resuelve un sistema de ecuaciones lineales triangular inferior Lx = b para x.
     *
     * @ N       El tamano de la matriz triangular inferior L y los vectores b y x.
     * @ L       La matriz triangular inferior de tamano N x N.
     * @ b       El vector de terminos independientes de tamano N.
     * @ x       El vector de solucion, que se actualizara con la solucion despues de la llamada.
     *
     * Esta funcion utiliza sustitucion hacia adelante para resolver el sistema de ecuaciones lineales.
     * La matriz L debe ser triangular inferior, y el vector x debe contener una estimacion inicial
     * de la solucion, que se actualizara con la solucion final.
     */

    double acumulador;
    for (int i = 0; i < N; i++) {
        acumulador = 0;
        for (int k = 0; k <= i; k++)
            acumulador += L[i][k] * x[k];
        x[i] = (b[i] - acumulador) / L[i][i];
    }
}

