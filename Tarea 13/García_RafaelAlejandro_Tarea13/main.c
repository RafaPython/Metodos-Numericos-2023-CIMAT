/*
Tarea 13 - Métodos Numéricos
Nombre: Rafael Alejandro García Ramírez
email: rafael.ramirez@cimat.mx
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "solvers.h"


double f1(double t, double y){
    return y;
}

int main(){


    printf("\n\n ---- METODOS CLASICOS DE SOLUCION DE EDO's ----\n\n");
    printf("Para el sistema y' = y\n");
    printf("1. Metodo de Euler (1 ecuacion)\n");
    printf("2. Metodo de Heun (1 ecuacion)\n");
    printf("3. Metodo de Taylor de segundo orden (1 ecuacion)\n");
    printf("4. Metodo de Runge - Kutta 4 (1 ecuacion)\n");
    printf("Para el problema de Lotka-Volterra\n");
    printf("5. Metodo de Euler (sistema de ecuaciones)\n");
    printf("6. Metodo de Heun (sistema de ecuaciones)\n");
    printf("7. Metodo de Taylor de segundo orden (sistema de ecuaciones)\n");
    printf("8. Metodo de Runge - Kutta 4 (sistema de ecuaciones)\n");
    printf("9. Salir\n");
    printf("Introduzca el numero del metodo que desea utilizar: ");
    int opcion;
    scanf("%d", &opcion);
    printf("\n");

    if (opcion < 1 || opcion > 9){
        printf("Opcion no valida\n");
        return 0;
    }

    if (opcion == 1){
        printf("Metodo de Euler (1 ecuacion)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, y0;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial: ");
        scanf("%lf", &y0);
        printf("\n");
        double *y_euler = euler(N, a, b, y0, f1);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("datos_euler.txt", "w");
        fprintf(fp, "t, euler\n");
        double copy = a;
        double h = (b-a)/N;
        for(int i = 0; i <= N; i++){
            fprintf(fp, "%f %f\n", copy, y_euler[i]);
            copy += h;
        }
        copy = a;
        fclose(fp);
        free(y_euler);
        printf("Los datos de la solucion se han guardado en el archivo datos_euler.txt\n");
    }
    if (opcion == 2){
        printf("Metodo de Heun (1 ecuacion)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, y0;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial: ");
        scanf("%lf", &y0);
        printf("\n");
        double *y_heun = heun(N, a, b, y0, f1);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("datos_heun.txt", "w");
        fprintf(fp, "t, heun\n");
        double copy = a;
        double h = (b-a)/N;
        for(int i = 0; i <= N; i++){
            fprintf(fp, "%f %f\n", copy, y_heun[i]);
            copy += h;
        }
        copy = a;
        fclose(fp);
        free(y_heun);
        printf("Los datos de la solucion se han guardado en el archivo datos_heun.txt\n");
    }
    if (opcion == 3){
        printf("Metodo de Taylor de segundo orden (1 ecuacion)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, y0;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial: ");
        scanf("%lf", &y0);
        printf("\n");
        double *y_taylor = taylor_segundo_orden(N, a, b, y0, f1);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("datos_taylor.txt", "w");
        fprintf(fp, "t, taylor\n");
        double copy = a;
        double h = (b-a)/N;
        for(int i = 0; i <= N; i++){
            fprintf(fp, "%f %f\n", copy, y_taylor[i]);
            copy += h;
        }
        copy = a;
        fclose(fp);
        free(y_taylor);
        printf("Los datos de la solucion se han guardado en el archivo datos_taylor.txt\n");
    }
    if (opcion == 4){
        printf("Metodo de Runge - Kutta 4 (1 ecuacion)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, y0;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial: ");
        scanf("%lf", &y0);
        printf("\n");
        double *y_RK4 = RK4(N, a, b, y0, f1);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("datos_RK4.txt", "w");
        fprintf(fp, "t, RK4\n");
        double copy = a;
        double h = (b-a)/N;
        for(int i = 0; i <= N; i++){
            fprintf(fp, "%f %f\n", copy, y_RK4[i]);
            copy += h;
        }
        copy = a;
        fclose(fp);
        free(y_RK4);
        printf("Los datos de la solucion se han guardado en el archivo datos_RK4.txt\n");
    }
    if (opcion == 5){
        printf("Metodo de Euler (sistema de ecuaciones)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, x0, y01;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial de x: ");
        scanf("%lf", &x0);
        printf("\n");
        printf("Introduzca el valor inicial de y: ");
        scanf("%lf", &y01);
        printf("\n");
        double y0[] = {x0, y01};
        int dim = 2;  
        double *resultados = euler_sistema(N, a, b, y0, dim, lotka_volterra);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("lotka_volterra_euler.txt", "w");
        fprintf(fp, "t x y\n");
        for(int i = 0; i <= N; i++) {
            fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
        }
        fclose(fp);
        free(resultados);
        printf("Los datos de la solucion se han guardado en el archivo lotka_volterra_euler.txt\n");
    }
    if (opcion == 6){
        printf("Metodo de Heun (sistema de ecuaciones)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, x0, y01;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial de x: ");
        scanf("%lf", &x0);
        printf("\n");
        printf("Introduzca el valor inicial de y: ");
        scanf("%lf", &y01);
        printf("\n");
        double y0[] = {x0, y01};
        int dim = 2;  
        double *resultados = heun_sistema(N, a, b, y0, dim, lotka_volterra);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("lotka_volterra_heun.txt", "w");
        fprintf(fp, "t x y\n");
        for(int i = 0; i <= N; i++) {
            fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
        }
        fclose(fp);
        free(resultados);
        printf("Los datos de la solucion se han guardado en el archivo lotka_volterra_heun.txt\n");
    }
    if (opcion == 7){
        printf("Metodo de Taylor de segundo orden (sistema de ecuaciones)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, x0, y01;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial de x: ");
        scanf("%lf", &x0);
        printf("\n");
        printf("Introduzca el valor inicial de y: ");
        scanf("%lf", &y01);
        printf("\n");
        double y0[] = {x0, y01};
        int dim = 2;  
        double *resultados = taylor_segundo_orden_sistema(N, a, b, y0, dim, lotka_volterra);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("lotka_volterra_taylor.txt", "w");
        fprintf(fp, "t x y\n");
        for(int i = 0; i <= N; i++) {
            fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
        }
        fclose(fp);
        free(resultados);
        printf("Los datos de la solucion se han guardado en el archivo lotka_volterra_taylor.txt\n");
    }
    if (opcion == 8){
        printf("Metodo de Runge - Kutta 4 (sistema de ecuaciones)\n");
        printf("Introduzca el numero de puntos: ");
        int N;
        scanf("%d", &N);
        printf("\n");
        printf("Introduzca el extremo inferior del intervalo: ");
        double a, b, x0, y01;
        scanf("%lf", &a);
        printf("\n");
        printf("Introduzca el extremo superior del intervalo: ");
        scanf("%lf", &b);
        printf("\n");
        printf("Introduzca el valor inicial de x: ");
        scanf("%lf", &x0);
        printf("\n");
        printf("Introduzca el valor inicial de y: ");
        scanf("%lf", &y01);
        printf("\n");
        double y0[] = {x0, y01};
        int dim = 2;  
        double *resultados = RK4_sistema(N, a, b, y0, dim, lotka_volterra);
        // guardamos los datos en un archivo
        FILE *fp;
        fp = fopen("lotka_volterra_RK4.txt", "w");
        fprintf(fp, "t x y\n");
        for(int i = 0; i <= N; i++) {
            fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
        }
        fclose(fp);
        free(resultados);
        printf("Los datos de la solucion se han guardado en el archivo lotka_volterra_RK4.txt\n");
    }
    if (opcion == 9){
        printf("Saliendo...\n");
        return 0;
    }

    return 0;
}
