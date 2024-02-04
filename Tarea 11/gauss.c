#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"


double cuadratura_gaussiana(double a, double b, double(*func)(double), int n){
    if (n < 1 || n > 4) {
        printf("\n---->Error: n debe ser un entero entre 1 y 4\n");
        return 0.0;
    }
    if (n == 1){
        return 2 * func((a + b) / 2)* (b - a) / 2;
    }
    else if (n == 2){
        double coeficientes[] = {1.0, 1.0};
        double raices[] = {0.5773502692, -0.5773502692};
        double resultado = 0.0;
        double cambio_limites;
        for(int i = 0; i < 2; i++) {
            cambio_limites = ((b-a)*raices[i] + (b+a))/2;
            resultado += coeficientes[i]*func(cambio_limites);
        }
        return resultado*(b - a) / 2.0;
    }
    else if (n == 3){
        double coeficientes[] = {0.555555555555556, 0.888888888888889, 0.555555555555556};
        double raices[] = {-0.774596669241483, 0.0, 0.774596669241483};
        double resultado = 0.0;
        double cambio_limites;
        for(int i = 0; i < 3; i++){
            cambio_limites = ((b-a)*raices[i] + (b+a))/2;
            resultado += coeficientes[i]*func(cambio_limites);
        }
        return resultado*(b - a) / 2.0;
    }
    else if (n == 4){
        double coeficientes[] = {0.652145154862546, 0.652145154862546, 0.347854845137454, 0.347854845137454};
        double raices[] = {-0.339981043584856, 0.339981043584856, -0.861136311594053, 0.861136311594053};
        double resultado = 0.0;
        double cambio_limites;
        for(int i = 0; i < 4; i++){
            cambio_limites = ((b-a)*raices[i] + (b+a))/2;
            resultado += coeficientes[i]*func(cambio_limites);
        }
        return resultado*(b - a) / 2.0;
    }
}


int main(){



    return 0;
}


