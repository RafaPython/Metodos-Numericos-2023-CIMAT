#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x){
    return 1.0/x;
}

double newton_cotes_cerrado(double(*func)(double), double a, double b, int n){
    if (n < 1 || n > 4) {
       printf("\n---->Error: n debe ser un entero entre 1 y 5\n");
       return 0.0;
    }
    double h = (b-a)/n;
    double resultado;
    switch (n) {
        case 1:
            resultado = 0.5*h*(func(a) + func(b));
            break;
        case 2:
            resultado = (h/3.0)*(func(a) + 4.0*func(a + h) + func(b));
            break;
        case 3:
            resultado = (3.0*h/8.0)*(func(a) + 3.0*func(a + h) + 3.0*func(a + 2.0*h) + func(b));
            break;
        case 4:
            resultado =(2.0*h/45.0)*(7.0*func(a) + 32.0*func(a + h) + 12.0*func(a + 2.0*h) + 32.0*func(a + 3.0*h) + 7.0*func(b));
            break;
    }
    return resultado;
}

double newton_cotes_abierto(double(*func)(double), double a, double b, int n){
    if (n < 0 || n > 3) {
       printf("\n---->Error: n debe ser un entero entre 0 y 3\n");
       return 0.0;
    }
    double h = (b-a)/(n+2.0);
    double resultado;
    switch (n) {
    case 0:
        resultado = 2*h*func(a + h);
        break;
    case 1:
        resultado = (3.0*h/2.0)*(func(a + h) + func(a + 2.0*h));
        break;
    case 2:
        resultado = (4.0*h/3.0)*(2.0*func(a + h) - func(a + 2.0*h) + 2.0*func(a + 3.0*h));
        break;
    case 3:
        resultado = (5.0*h/24.0)*(11.0*func(a + h) + func(a + 2.0*h) + func(a + 3.0*h) + 11.0*func(a + 4.0*h));
        break;
    }
    return resultado;
}
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
    return 0.0;
}

int main(){

    printf("Problema 1\n");

    double a = 1.0;
    double b = 2.0;

    double real = log(2.0);
    double nc, na, cg;

    for (int i = 2; i <= 3; i++){   
        printf("\n");
        printf("n = %d\n", i);
        nc = newton_cotes_cerrado(f, a, b, i);
        na = newton_cotes_abierto(f, a, b, i);
        cg = cuadratura_gaussiana(a, b, f, i);
        printf("Newton-Cotes cerrado: %.15lf\n", nc);
        printf("El error absoluto es: %.15lf\n", fabs(real - nc));
        printf("Newton-Cotes abierto: %.15lf\n", na);
        printf("El error absoluto es: %.15lf\n", fabs(real - na));
        printf("Cuadratura Gaussiana: %.15lf\n", cg);
        printf("El error absoluto es: %.15lf\n", fabs(real - cg));
    
        printf("\n");

    }


    return 0;
}