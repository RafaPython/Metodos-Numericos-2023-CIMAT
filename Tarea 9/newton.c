#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

double interpolador_newton(int N, double *x, double *y, double punto, double **tabla);

int main(){
    int N = 4;
    double *x = (double*)malloc(N*sizeof(double));
    double *y = (double*)malloc(N*sizeof(double));
    double **tabla = generar_matriz(N);

    x[0] = -1.0; x[1] = 0.0; x[2] = 1.0; x[3] = 2.0;
    y[0] = 3.0; y[1] = -4.0; y[2] = 5.0; y[3] = -6.0;

    double punto = 0.0;

    printf("El punto es %lf\n\n", punto);

    // imprimitmos los datos
    printf("Los datos ingresados son:\n");
    for(int i = 0; i < N; i++){
        printf("x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
    }
    printf("\n");

    double P = interpolador_newton(N, x, y, punto, tabla);

    printf("El valor de la interpolacion es: %.12lf\n", P);

    // tabla de diferencias divididas
    printf("\n");
    printf("La tabla de diferencias divididas es:\n");
    imprimir_matriz(N, tabla);

    free(x);
    free(y);
    liberar_matriz(N, tabla);



    return 0;
}

double interpolador_newton(int N, double *x, double *y, double punto, double **tabla){
    // hacemos la tabla de diferencias divididas cero 
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) tabla[i][j] = 0.0;


    // inicializamos la tabla de diferencias divididas
    for(int i = 0; i < N; i++) tabla[i][0] = y[i];

    // calculamos las diferencias divididas
    for(int i = 1; i < N; i++){
        for(int j = 1; j < i + 1; j++){
            tabla[i][j] = (tabla[i][j-1] - tabla[i-1][j-1])/(x[i] - x[i-j]);
        }
    }

    // calculamos el polinomio interpolador
    double L, P = 0.0;
    for(int i = 1; i < N; i++){
        L = 1;
        for(int j = 0; j < N; j++){
            if(j <= i - 1) L *= (punto - x[j]);
        }
        P += tabla[i][i] * L;
    }

    return P + tabla[0][0];
}