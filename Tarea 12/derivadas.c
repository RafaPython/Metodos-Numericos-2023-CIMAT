#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"

/*
Programa que realiza las siguiente tipos de derivadas

- Derivada hacia adelante
- Derivada hacia atras
- Derivada centrada
- Derivada de tres puntos
- Derivada de cinco puntos
- Derivada de segundo orden

*/
/// formulas generales 
// double derivada_hacia_adelante(double (*f)(double), double x, double h){
//     return (f(x+h)-f(x))/h;
// }
// double derivada_hacia_atras(double (*f)(double), double x, double h){
//     return (f(x)-f(x-h))/h;
// }
// double derivada_centrada(double (*f)(double), double x, double h){
//     return (f(x+h)-f(x-h))/(2*h);
// }

// /// /////////

// double derivada_tres_puntos_endpoint(double (*f)(double), double x, double h){
//     return (1/(2*h))*(-3*f(x)+4*f(x+h)-f(x+2*h));
// }

// double derivada_tres_puntos_midpoint(double (*f)(double), double x, double h){
//     return (1/(2*h))*(f(x+h)-f(x-h));
// }

// double derivada_cinco_puntos_midpoint(double (*f)(double), double x, double h){
//     return (1/(12*h))*(f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h));
// }

// double derivada_cinco_puntos_endpoint(double (*f)(double), double x, double h){
//     return (1/(12*h))*(-25*f(x) + 48*f(x+h) - 36*f(x+2*h) + 16*f(x+3*h) - 3*f(x+4*h));
// }
// double segunda_derivada(double (*f)(double), double x, double h){
//     return (1/(h*h))*(f(x-h)-2*f(x)+f(x+h));
// }
// // ahora lo hacemos para datos en lugar de funciones

// double derivada_hacia_adelante_datos(double f1, double f2, double h){
//     // f1 es el dato en f(x)
//     // f2 es el dato en f(x+h)
//     return (f2-f1)/h;
// }
// double derivada_hacia_atras_datos(double f1, double f2, double h){
//     // f1 es el dato en f(x-h)
//     // f2 es el dato en f(x)
//     return (f2-f1)/h;
// }

// double derivada_centrada_datos(double f1, double f2, double h){
//     // f1 es el dato en f(x-h)
//     // f2 es el dato en f(x+h)
//     return (f2-f1)/(2*h);
// }

// double derivada_tres_puntos_endpoint_datos(double f1, double f2, double f3, double h){
//     // f1 es el dato en f(x)
//     // f2 es el dato en f(x+h)
//     // f3 es el dato en f(x+2h)
//     return (1/(2*h))*(-3*f1+4*f2-f3);
// }
// double derivada_tres_puntos_midpoint_datos(double f1, double f2, double h){
//     // f1 es el dato en f(x+h)
//     // f2 es el dato en f(x-h)
//     return (1/(2*h))*(f1-f2);
// }
// double derivada_cinco_puntos_midpoint_datos( double f1, double f2, double f3, double f4, double h){
//     // f1 es el dato en f(x-2h)
//     // f2 es el dato en f(x-h)
//     // f3 es el dato en f(x+h)
//     // f4 es el dato en f(x+2h)
//     return (1/12*h)*(f1-8*f2+8*f3-f4);
// }
// double derivada_cinco_puntos_endpoint_datos(double f1, double f2, double f3, double f4, double f5, double h){
//     // f1 es el dato en f(x)
//     // f2 es el dato en f(x+h)
//     // f3 es el dato en f(x+2h)
//     // f4 es el dato en f(x+3h)
//     // f5 es el dato en f(x+4h)
//     return (1/(12*h))*(-25*f1+48*f2-36*f3+16*f4-3*f5);
// }
// double segunda_derivada_midpoint_datos(double f1, double f2, double f3, double h){
//     // f1 es el dato en f(x-h)
//     // f2 es el dato en f(x)
//     // f3 es el dato en f(x+h)
//     return (1/(h*h))*(f1-2*f2+f3);

// }

double funcion_derivada(double x){
    return sin(x) + (3*exp(x)*(x+1));
}
double funcion_segunda_derivada(double x){
    return cos(x) + (3*exp(x)*(x+2));
}

int main(){
    // double x[] = {1.20, 1.29, 1.30, 1.31, 1.40};
    double f[] = {11.59006, 13.78176, 14.04276, 14.30741, 16.86187};

    double h1 = 0.1;
    double h2 = 0.01;
    double punto = 1.3;
    double real = funcion_derivada(punto);
    double real2 = funcion_segunda_derivada(punto);
    printf("Para h = 0.01\n");
    // El valor de la derivada hacia adelante en el punto 1.3 correspendo a f[3] - f[2] / h2
    printf("El valor de la derivada hacia adelante para x = 1.3 es: %lf\n", derivada_hacia_adelante_datos(f[2], f[3], h2));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error absoluto de la derivada hacia adelante para x = 1.3 es: %lf\n", fabs(derivada_hacia_adelante_datos(f[2], f[3], h2) - real));
    printf("\n");
    // el valor de la derivada hacia atrás será f[1] - f[2] / h2
    printf("El valor de la derivada hacia atrás para x = 1.3 es: %lf\n", derivada_hacia_atras_datos(f[1], f[2], h2));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error absoluto de la derivada hacia atrás para x = 1.3 es: %lf\n", fabs(derivada_hacia_atras_datos(f[1], f[2], h2) - real));
    printf("\n");
    // el valor de derivada centrada será f[3] - f[1] / 2h2
    printf("El valor de la derivada centrada para x = 1.3 es: %lf\n", derivada_centrada_datos(f[1], f[3], h2));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error absoluto de la derivada centrada para x = 1.3 es: %lf\n", fabs(derivada_centrada_datos(f[1], f[3], h2) - real));
    printf("\n");
    // el valor de la derivada de tres puntos midpoint será f[3] - f[1] / 2h2
    printf("El valor de la derivada de tres puntos midpoint para x = 1.3 es: %lf\n", derivada_tres_puntos_midpoint_datos(f[3], f[1], h2));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error absoluto de la derivada de tres puntos midpoint para x = 1.3 es: %lf\n", fabs(derivada_tres_puntos_midpoint_datos(f[3], f[1], h2) - real));
    printf("\n");

    //  los demás no se pueden dado que no se nos otorgan los puntos x para poder calcularlos

    // el valor de la segunda derivada midpoint es f[1] - 2*f[2] + f[3] / h2^2
    printf("El valor de la segunda derivada midpoint para x = 1.3 es: %lf\n", segunda_derivada_midpoint_datos(f[1], f[2], f[3], h2));
    // printf("El valor real de la segunda derivada para x = 1.3 es: %lf\n", funcion_segunda_derivada(punto));
    printf("El error absoluto de la segunda derivada midpoint para x = 1.3 es: %lf\n", fabs(segunda_derivada_midpoint_datos(f[1], f[2], f[3], h2) - real2));

    printf("\nPara h = 0.1\n");
    // El valor de la derivada hacia adelante en el punto 1.3 correspendo a f[4] - f[2] / h2
    printf("El valor de la derivada hacia adelante para x = 1.3 es: %lf\n", derivada_hacia_adelante_datos(f[2], f[4], h1));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error absoluto de la derivada hacia adelante para x = 1.3 es: %lf\n", fabs(derivada_hacia_adelante_datos(f[2], f[4], h1) - real));
    printf("\n");
    // el valor de la derivada hacia atrás será f[0] - f[2] / h2
    printf("El valor de la derivada hacia atrás para x = 1.3 es: %lf\n", derivada_hacia_atras_datos(f[0], f[2], h1));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error absoluto de la derivada hacia atrás para x = 1.3 es: %lf\n", fabs(derivada_hacia_atras_datos(f[0], f[2], h1) - real));
    printf("\n");
    // el valor de derivada centrada será f[4] - f[0] / 2h2
    printf("El valor de la derivada centrada para x = 1.3 es: %lf\n", derivada_centrada_datos(f[0], f[4], h1));
    // printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("El error para la derivada centrada para x = 1.3 es: %lf\n", fabs(derivada_centrada_datos(f[0], f[4], h1) - real));
    printf("\n");
    // el valor de la derivada de tres puntos midpoint será f[3] - f[1] / 2h2
    printf("El valor de la derivada de tres puntos midpoint para x = 1.3 es: %lf\n", derivada_tres_puntos_midpoint_datos(f[4], f[0], h1));
    printf("El valor real de la derivada para x = 1.3 es: %lf\n", funcion_derivada(punto));
    printf("\n");

    printf("El valor de la segunda derivada midpoint para x = 1.3 es: %lf\n", segunda_derivada_midpoint_datos(f[0], f[2], f[4], h1));
    printf("El valor real de la segunda derivada para x = 1.3 es: %lf\n", funcion_segunda_derivada(punto));
    printf("\n");


    return 0;
}