#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solvers.h"
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
        printf("t[%d] = %f, y[%d] = %f\n", i, x,i, y[i]);
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
        printf("t[%d] = %f, y[%d] = %f\n", i, x ,i, y[i]);
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
        printf("x[%d] = %f, y[%d] = %f\n", i, x ,i, y[i]);
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
        printf("t[%d] = %f, y[%d] = %f\n", i, x ,i, y[i]);
        x += h;
    }
    return y;
}
    


double fe(double x, double y){
    return sqrt(1 + (x*x*x));
}
double f1(double x){
    return sqrt(1 + (x*x*x));
}

int main(){

    double a = 0;
    double b = 2;
    double I = cuadratura_gaussiana(a, b, f1, 3);
    printf("I = %.20f\n", I);
    double y0 = 0.0;
    int N = 10;
    double *y_euler = euler(N, a, b, y0, fe);
    double *y_heun = heun(N, a, b, y0, fe);
    double *y_taylor = taylor_segundo_orden(N, a, b, y0, fe);
    double *y_RK4 = RK4(N, a, b, y0, fe);

    // comparamos los resultados con la solucion analitica
    printf("\nError relativo de la solucion analitica\n");
    printf("Error con respecto al metodo de Euler: %f\n", fabs(y_euler[N] - I));
    printf("Error con respecto al metodo de Heun: %f\n", fabs(y_heun[N] - I));
    printf("Error con respecto al metodo de Taylor de segundo orden: %f\n", fabs(y_taylor[N] - I));
    printf("Error con respecto al metodo de Runge-Kutta de orden 4: %f\n", fabs(y_RK4[N] - I));


    // ahora calculamos el valor de la integral con  h = 0.5 para el método de taylor de segundo orden
    N = 4;
    printf("\nPara h = 0.5\n");
    double *y_taylor2 = taylor_segundo_orden(N, a, b, y0, fe);
    printf("I = %.20f\n", y_taylor2[N]);

    // guardamos los datos en un archivo
    FILE *fp;
    fp = fopen("datos_integral.txt", "w");
    fprintf(fp, "t, euler, heun, taylor, RK4\n");
    double copy = a;
    double h = (b-a)/N;
    for(int i = 0; i <= N; i++){
        fprintf(fp, "%f %f %f %f %f\n", copy, y_euler[i], y_heun[i], y_taylor[i], y_RK4[i]);
        copy += h;
    }
    copy = a;
    fclose(fp);


    free(y_euler);
    free(y_heun);
    free(y_taylor);
    free(y_RK4);

    return 0;
}