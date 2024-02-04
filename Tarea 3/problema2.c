/*
Problema 2 Tarea 3 Métodos Numéricos 

Autor: Rafael Alejandro García Ramírez
email: rafael.ramirez@cimat.mx

*/


#include <stdio.h>

void Ux_b(int N, double U[N][N], double b[N] ,double x[N]);

int main(){
    // Ejemplo
    const int N = 3;
    double U[3][3] = {{3,1,3}, {0,4,3}, {0,0,5}};
    double b[3] = {2,5,8};

    double x[N];
    Ux_b(N, U, b, x);
    for (int i = 0; i < N; i++) printf("x[%d] = %f\n", i, x[i]);
    return 0;
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
        for (int k =  N - 1; k > i; k--) acumulador += U[i][k]*x[k]; // N - 1 or N
        x[i] = (b[i] - acumulador) / U[i][i];
    }  
}

