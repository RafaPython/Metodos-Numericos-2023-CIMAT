#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"

/*
Programa que resuelve el problema de minimos cuadrados para un polinomio de grado n que interpola n+1 puntos dados.
P(x) = a0 + a1*x + a2*x^2 + ... + an*x^n

*/
void interpolacion_cuadrados_polinomio(int N, double *x, double *y, double *w, double regulador);

int main(){

    // Definimos los puntos a interpolar
    int N = 5;
    double *x = (double *)malloc(N*sizeof(double));
    double *y = (double *)malloc(N*sizeof(double));
    double *w = (double *)malloc(N*sizeof(double));
    double regulador = 0.0;

    x[0] = 0.0; y[0] = 1.0;
    x[1] = 1.0; y[1] = 2.0;
    x[2] = 2.0; y[2] = 3.0;
    x[3] = 2.4; y[3] = 4.2;
    x[4] = 4.0; y[4] = 10.0;

    // Calculamos los pesos
    interpolacion_cuadrados_polinomio(N, x, y, w, regulador);

    // Imprimimos los resultados
    printf("\n\n");
    for(int i = 0; i < N; i++) printf("w[%d] = %lf\n", i, w[i]);
    printf("\n\n");

    // definimos el punto a evaluar
    double x0 = 1.5;
    double y0 = 0.0;
    for(int i = 0; i < N; i++) y0 += w[i]*pow(x0, i);
    printf("P(%lf) = %lf\n", x0, y0);

    // liberamos la memoria
    free(x);
    free(y);
    free(w);

    return 0;
}

void interpolacion_cuadrados_polinomio(int N, double *x, double *y, double *w, double regulador){

    /*
        @ N: numero de puntos
        @ x: vector de abscisas
        @ y: vector de ordenadas
        @ w: vector de pesos (vacío para escribir los pesos calculados)

    */

    // Creamos la matriz de diseño phi  phi_k = sum phi^T * phi 

    double **phi = generar_matriz(N);
    for (int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            phi[i][j] = pow(x[i], j);
    
    // realizamos el producto phi^T * y
    double **phiT = generar_matriz(N);
    trasponer_matriz(N, phi, phiT);

    // realizamos el producto phi^T * phi
    double **phiTphi = generar_matriz(N);
    multiplicar_matrices(N, phiT, phi, phiTphi);
    for(int i = 0; i < N; i++) phiTphi[i][i] += regulador;

    // multiplicamos phi^T * y
    double *phiTy = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, phiTy);
    


    double *d = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, d);

    // Resolvemos el sistema phiTphi * w = d
    // rellenamos w con valores aleatorios
    double **L = generar_matriz(N);
    double **U = generar_matriz(N);
    doolittle(N, phiTphi, L, U);


    double *dummy = (double *)malloc(N*sizeof(double));
    // resolvemos el sistema L * dummy = d
    Lx_b(N, L, d, dummy);
    // resolvemos el sistema U * w = dummy
    Ux_b(N, U, dummy, w);

    // liberamos la memoria
    liberar_matriz(N, phi);
    liberar_matriz(N, phiT);
    liberar_matriz(N, phiTphi);
    free(phiTy);
    free(d);
    free(dummy);
    liberar_matriz(N, L);
    liberar_matriz(N, U);
}

void interpolacion_cuadrados_coseno(int N, double *x, double *y, double *w, double regulador){

    /*
        @ N: numero de puntos
        @ x: vector de abscisas
        @ y: vector de ordenadas
        @ w: vector de pesos (vacío para escribir los pesos calculados)

    */

    // Creamos la matriz de diseño phi  phi_k = sum phi^T * phi 
    double PI = 3.14159265358979323846;
    double **phi = generar_matriz(N);
    for (int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            phi[i][j] = cos(j*x[i]*PI/6.0);
    
    // realizamos el producto phi^T * y
    double **phiT = generar_matriz(N);
    trasponer_matriz(N, phi, phiT);

    // realizamos el producto phi^T * phi
    double **phiTphi = generar_matriz(N);
    multiplicar_matrices(N, phiT, phi, phiTphi);
    for(int i = 0; i < N; i++) phiTphi[i][i] += regulador;

    // multiplicamos phi^T * y
    double *phiTy = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, phiTy);
    


    double *d = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, d);

    // Resolvemos el sistema phiTphi * w = d
    // rellenamos w con valores aleatorios
    double **L = generar_matriz(N);
    double **U = generar_matriz(N);
    doolittle(N, phiTphi, L, U);


    double *dummy = (double *)malloc(N*sizeof(double));
    // resolvemos el sistema L * dummy = d
    Lx_b(N, L, d, dummy);
    // resolvemos el sistema U * w = dummy
    Ux_b(N, U, dummy, w);

    // liberamos la memoria
    liberar_matriz(N, phi);
    liberar_matriz(N, phiT);
    liberar_matriz(N, phiTphi);
    free(phiTy);
    free(d);
    free(dummy);
    liberar_matriz(N, L);
    liberar_matriz(N, U);
}

void interpolacion_cuadrados_exponencial(int N, double *x, double *y, double *w, double regulador){ 
        /*
            @ N: numero de puntos
            @ x: vector de abscisas
            @ y: vector de ordenadas
            @ w: vector de pesos (vacío para escribir los pesos calculados)
    
        */

        // Creamos la matriz de diseño phi  phi_k = sum phi^T * phi
        double **phi = generar_matriz(N);
        double r;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++){ 
                r = (x[i] - x[j])*(x[i] - x[j]);
                phi[i][j] = exp(-r);
            }
                
        
        // realizamos el producto phi^T * y
        double **phiT = generar_matriz(N);
        trasponer_matriz(N, phi, phiT);

        // realizamos el producto phi^T * phi
        double **phiTphi = generar_matriz(N);
        multiplicar_matrices(N, phiT, phi, phiTphi);
        for (int i = 0; i < N; i++)
            phiTphi[i][i] += regulador;

        // multiplicamos phi^T * y
        double *phiTy = (double *)malloc(N * sizeof(double));
        multiplicar_matriz_vector(N, phiT, y, phiTy);

        double *d = (double *)malloc(N * sizeof(double));
        multiplicar_matriz_vector(N, phiT, y, d);

        // Resolvemos el sistema phiTphi * w = d
        // rellenamos w con valores aleatorios
        double **L = generar_matriz(N);
        double **U = generar_matriz(N);
        doolittle(N, phiTphi, L, U);

        double *dummy = (double *)malloc(N * sizeof(double));
        // resolvemos el sistema L * dummy = d
        Lx_b(N, L, d, dummy);
        // resolvemos el sistema U * w = dummy
        Ux_b(N, U, dummy, w);

        // liberamos la memoria
        liberar_matriz(N, phi);
        liberar_matriz(N, phiT);
        liberar_matriz(N, phiTphi);
        free(phiTy);
        free(d);
        free(dummy);
        liberar_matriz(N, L);
        liberar_matriz(N, U);
}