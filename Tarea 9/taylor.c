#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

double interpolacion_taylor(int N, double *derivadas, double x, double x0);

double error_absoluto(double valor_real, double valor_aproximado){
    return fabs(valor_real - valor_aproximado);
}

int main(){
    // function f(x) = 1 / (1 + x^2)

    int N = 4;
    int n[] = {1,3,5,10};
    

    double x0 = 0;
    double x[] = {0.5,1.0, 1.5, 2.0};

    double valor_real[] = {exp(0.5), exp(1.0), exp(1.5), exp(2.0)};

    // double derivadas[] = {1, 1, 1, 1}; // para el caso e^x


    // imprimmos los datos 
    // printf("Los datos ingresados son:\n");
    // for (int i = 0; i < N; i++){
    //     printf("derivadas[%d] = %lf\n", i, derivadas[i]);
    // }
    // printf("\n");
    // printf("x0 = %lf\n", x0);

    // iteramos sobre los datos de x
    for(int i = 0; i < N; i++){
        printf("x = %lf\n", x[i]);
        double *derivadas = (double*)malloc(n[i]*sizeof(double));
        for(int k = 0; k < n[i]; k++) derivadas[k] = 1;

        for(int j = 0; j < 4; j++){
            printf("\nEl orden de la interpolacion es: %d\n", n[j]);
            double P = interpolacion_taylor(n[j], derivadas, x[i], x0);
            printf("El valor de la interpolacion es: %.12lf\n", P);
            printf("El error absoluto es: %.12lf\n", error_absoluto(valor_real[i], P));
        }
        printf("\n");
        free(derivadas);
    }



    return 0;
}

// double interpolacion_taylor(int N, double *derivadas, double x, double x0){
//     double diferencia = x - x0;
//     double diferencia_previa = 1;
//     double factorial = 1;

//     double P = derivadas[0];

//     for(int i = 1; i < N ; i++){
//         diferencia_previa *= diferencia;
//         factorial *= i ;
//         P += (derivadas[i] * diferencia_previa) / factorial;  
//     }

//     return P;
// }