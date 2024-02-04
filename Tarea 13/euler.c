#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"

/*
Programa para resolver una ecuacion diferencial ordinaria mediante el método de Euler, el método de Heun y el método de Runge-Kutta de orden 4.
*/

double *euler(int N, double a, double b, double y0, double (*f)(double, double)){
    double h = (b - a) / N;
    double *y = malloc((N+1)*sizeof(double));
    double x = a;
    y[0] = y0;
    for(int i = 1; i <= N; i++){
        y[i] = y[i-1] + h*f(a, y[i-1]);
        a += h;
    }
    // imprimimos los datos
    printf("\nSolucion por Euler\n");
    for(int i = 0; i <= N; i++){
        printf("t[%d] = %f, y[%d] = %f ", i, x,i, y[i]);
        printf("    e_r = %f\n", fabs(y[i] - exp(x)));
        x += h;
    }
    return y;
}


double *heun(int N, double a, double b, double y0, double (*f)(double, double)){
    double h = (b - a) / N;
    double *y = malloc((N+1)*sizeof(double));
    double x = a;
    double yDummy, a_anterior = a;
    y[0] = y0;
    for(int i = 1; i <= N; i++){
        a += h;
        yDummy = y[i-1] + h*f(a_anterior, y[i-1]);
        y[i] = y[i-1] + (h/2)*(f(a_anterior, y[i-1]) + f(a, yDummy));
        a_anterior = a;
    }
    // // imprimimos los datos
    printf("\nSolucion por Heun\n");
    for(int i = 0; i <= N; i++){
        printf("t[%d] = %f, y[%d] = %f", i, x ,i, y[i]);
        printf("   e_r =  %f\n", fabs(y[i] - exp(x)));
        x += h;
    }
    return y;
}


double derivada_parcial_x(double (*f)(double, double), double x, double y, double h, int variable){
    if(variable == 0){
        return (f(x+h, y) - f(x-h, y))/(2*h);
    }
    if (variable == 1){
        return (f(x, y+h) - f(x, y-h))/(2*h);
    }
    else{
        printf("Error en la variable\n");
        return 0;
    }
}

double *taylor_segundo_orden(int N, double a, double b, double y0, double (*f)(double, double)){
    double h = (b - a) / N;
    double *y = malloc((N+1)*sizeof(double));
    double derivada_x, derivada_y;
    double h_2 = (h*h) /2;
    y[0] = y0;
    double x = a;
    for(int i = 1; i <= N; i++){
        derivada_x = derivada_parcial_x(f, a, y[i-1], h, 0);
        derivada_y = derivada_parcial_x(f, a, y[i-1], h, 1);
        y[i] = y[i-1] + h*f(a, y[i-1]) + h_2*(derivada_x + derivada_y*f(a, y[i-1]));
        a += h;
    }
    // // imprimimos los datos
    printf("\nSolucion por Taylor de segundo orden\n");
    for(int i = 0; i <= N; i++){
        printf("t[%d] = %f, y[%d] = %f", i, x ,i, y[i]);
        printf("    e_r = %f\n", fabs(y[i] - exp(x)));
        x += h;
    }
    return y;
}

double *RK4(int N, double a, double b, double y0, double (*f)(double, double)) {
    double h = (b - a) / N;
    double *y = malloc((N+1)*sizeof(double));
    double x = a;
    double k1, k2, k3, k4;
    y[0] = y0;
    for(int i = 1; i <= N; i++){
        k1 = h*f(a, y[i-1]);
        k2 = h*f(a + h/2, y[i-1] + k1/2);
        k3 = h*f(a + h/2, y[i-1] + k2/2);
        k4 = h*f(a + h, y[i-1] + k3);
        y[i] = y[i-1] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
        a += h;
    }
    // // imprimimos los datos
    printf("\nSolucion por Runge-Kutta de orden 4\n");
    for(int i = 0; i <= N; i++){
        printf("t[%d] = %f, y[%d] = %f", i, x ,i, y[i]);
        printf("    e_r = %f\n", fabs(y[i] - exp(x)));
        x += h;
    }
    return y;
}
    


double f1(double t, double y){
    return y;

}
int main(){
    int N = 99;
    double a = 0.0;
    double b = 4.0;
    double y0 = 1.0;
    double *y = euler(N, a, b, y0, f1);

    // guardamos los datos dentro de un archivo 
    FILE *fp;
    fp = fopen("datos_euler.txt", "w");
    double copy = a;
    double h = (b-a)/N;
    for(int i = 0; i <= N; i++){
        fprintf(fp, "%f %f\n", copy, y[i]);
        copy += h;
    }
    copy = a;
    free(y);
    double *y2 = heun(N, a, b, y0, f1);
    fp = fopen("datos_heun.txt", "w");
    for(int i = 0; i <= N; i++){
        fprintf(fp, "%f %f\n", copy, y2[i]);
        copy += h;
    }
    copy = a;
    free(y2);
    double *y3 = taylor_segundo_orden(N, a, b, y0, f1);
    fp = fopen("datos_taylor.txt", "w");
    for(int i = 0; i <= N; i++){
        fprintf(fp, "%f %f\n", copy, y3[i]);
        copy += h;
    }
    copy = a;
    free(y3);
    double *y4 = RK4(N, a, b, y0, f1);
    fp = fopen("datos_RK4.txt", "w");
    for(int i = 0; i <= N; i++){
        fprintf(fp, "%f %f\n", copy, y4[i]);
        copy += h;
    }
    copy = a;
    free(y4);
    // creamos un archivo con los datos de la solucion analitica
    fp = fopen("datos_analitica.txt", "w");
    for(int i = 0; i <= N; i++){
        fprintf(fp, "%f %f\n", copy, exp(copy));
        copy += h;
    }
    return 0;
}