#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"


// Función para calcular el polinomio de Hermite
// double hermite(int N, double *x, double *y, double *derivadas, double punto) {
//     // abrimos espacio en la memoria 
//     double **Q = generar_matriz(N);

//     // inicializamos la matriz en ceros
//     for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) Q[i][j] = 0.0;

//     // llenamos la matriz con los valores de y y las derivadas
//     for(int i = 0; i < N; i++){
//         Q[i][0] = y[i];
//         Q[i][1] = derivadas[i];
//         if(i != 0) Q[i][1] = (Q[i][0] - Q[i - 1][0]) / (x[i] - x[i - 1]);
//     }

//     // calculamos los coeficientes
//     for(int i = 2; i < N; i++)
//          for (int j = 2; j <= i; j++)
//              Q[i][j] = (Q[i][j - 1] - Q[i - 1][j - 1]) / (x[i] - x[i - j]);
    
//     // calculamos el polinomio
//     double producto, P = 0.0;
//     for(int i = 0; i < N; i++){
//         producto = Q[i][i];
//         for(int j = 0; j < i; j++) producto *= (punto - x[j]);
//         P += producto;
//     }

//     // liberamos la memoria 
//     liberar_matriz(N, Q);
    
//     return P;
// }

int main() {

    int N = 4;
    double x[] ={0.30, 0.32, 0.33, 0.35};
    double y[] = {0.29552, 0.31457, 0.32404, 0.34290};
    double derivadas[] = {0.95534, 0.94924,0.94604, 0.93937};

    double punto = 0.34;

    // imprimir los puntos
    printf("Los puntos son: \n");
    for(int i = 0; i < N; i++) printf("(%lf, %lf)\n", x[i], y[i]);

    // imprimimos las derivadas
    printf("Las derivadas son: \n");
    for(int i = 0; i < N; i++) printf("f'(%lf) = %lf\n", x[i], derivadas[i]);

    // realizamos la interpolación de Hermite

    double P = hermite(N, x, y, derivadas, punto);

    printf("El valor de la interpolación de Hermite es: %lf\n", P);
    printf("El valor real es: %lf\n", sin(punto));
    double dummy = fabs(P - sin(punto));
    printf("El error absoluto es: %.20lf\n", dummy);

    printf("la condición con respecto a un punto menos es %s", dummy < 0.00000224119251895916 ? "verdadera" : "falsa");



    
    return 0;
}