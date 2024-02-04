#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// #include "solvers.h"

double interpolador_lagrange(int N, double *x, double *y, double punto);

int main(){

    int N = 5;

    double *x = (double*)malloc(N*sizeof(double));
    double *y = (double*)malloc(N*sizeof(double));

    // x[0] = 2; x[1] = 11.0/4; x[2] = 4;
    // y[0] = 1.0/2; y[1] = 4.0/11; y[2] = 1.0/4;

    x[0] = 1; x[1] = 1.3; x[2] = 1.6; x[3] = 1.9; x[4] = 2.2;
    y[0] = 0.1411; y[1] = -0.6878; y[2] = -0.9962; y[3] = -0.5507; y[4] = 0.3115;

    double punto = 1.5;

    printf("El punto es %lf\n\n", punto);

    // imprimimos los datos 
    printf("Los datos ingresados son:\n");
    for (int i = 0; i < N; i++){
        printf("x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
    }
    printf("\n");

    double P = interpolador_lagrange(N, x, y, punto);


    printf("El valor de la interpolacion es: %.12lf\n", P);

    free(x);
    free(y);
    
    return 0;
}

double interpolador_lagrange(int N, double *x, double *y, double punto){

    double L, P = 0.0;
    for (int i = 0; i < N; i++){
        L = 1;
        for (int j = 0; j < N; j++){
            if (i != j)
                 L *= (punto - x[j])/(x[i] - x[j]);
        }
        P += y[i] * L;
    }
    return P;
}



