#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

double interpolacion_neville(int N, double *x, double *y, double **Q, double punto);

int main(){
    int N = 5;
    double x[] = {1.0, 1.3, 1.6, 1.9, 2.2};
    double y[] = {0.7651977, 0.6200860, 0.4554022, 0.28148186, 0.1103626};
    double punto = 1.5;

    double **Q = generar_matriz(N);

    double P = interpolacion_neville(N, x, y, Q, punto);

    // imprimimos los datos 
    printf("Los datos ingresados son:\n");
    for (int i = 0; i < N; i++){
        printf("x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
    }
    printf("\n");
    printf("El punto es %lf\n\n", punto);

    printf("El valor de la interpolacion es: %.12lf\n", P);

    // tabla de diferencias divididas
    printf("\n");
    printf("La tabla de diferencias divididas es:\n");
    imprimir_matriz(N, Q);
    
    liberar_matriz(N, Q);

    return 0; 
}


double interpolacion_neville(int N, double *x, double *y, double **Q, double punto){

    // hacemos la tabla de diferencias divididas cero 
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) Q[i][j] = 0.0;

    // asignamos los valores de la tabla
    for(int i = 0; i < N; i++) Q[i][0] = y[i];

    // calculamos los valores de la tabla
    for(int i = 1; i < N; i++){
        for(int j = 1; j < i + 1; j++){
            Q[i][j] = ((punto - x[i-j])*Q[i][j-1] - (punto - x[i])*Q[i-1][j-1])/(x[i] - x[i-j]);
        }
    }

    return Q[N-1][N-1];
}
