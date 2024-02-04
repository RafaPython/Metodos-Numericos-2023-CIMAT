#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

/*
Nombre : Rafael Alejandro García Ramírez 
email : rafael.ramirez@cimat.mx

NOTA: realmente no sabía si debía de poner una opción de lectura de datos, si después 
se requiere es solamente agregar una parte para lectura de archivos con un 
formato de datos específico.

*/


int main(){

    printf("\n--- METODOS DE INTERPOLACION --- \n\n");

    int opcion = 0;
    printf("Seleccione el metodo de interpolacion:\n");
    printf("1. Taylor\n");
    printf("2. Lagrange\n");
    printf("3. Neville\n");
    printf("4. Newton\n");
    printf("5. Salir\n");
    printf("Opcion: ");
    scanf("%d", &opcion);
    
    if (opcion == 5) {
        printf("\n\nSaliendo...\n\n");
        return 0;
    }
    if (opcion < 1 || opcion > 5){
        printf("\n\nOpcion no valida\n\n");
        return 0;
    }

    int N;
    printf("\n\nIngrese el numero de datos con los que se desea calcular: ");
    scanf("%d", &N);

    double *x = (double*)malloc(N*sizeof(double));
    double *y = (double*)malloc(N*sizeof(double));
    double **tabla = generar_matriz(N);

    double x0, P, punto;

    switch(opcion){
        case 1: 
            // double interpolacion_taylor(int N, double *derivadas, double x, double x0);
            printf("\n\n--- METODO DE TAYLOR ---\n\n");

            printf("Ingrese el punto en el que se desea calcular: ");
            scanf("%lf", &punto);
            
            printf("\nIngrese el punto x0 sobre el cual se evalúan las derivadas: ");
            scanf("%lf", &x0);
            printf("\n");

            for(int i = 0; i < N; i++){
                printf("\nIngrese el valor de la derivada %d: ", i);
                scanf("%lf", &x[i]);
            }

            printf("\n\n");
            
            // imprimimos los datos
            printf("Los datos ingresados son:\n");
            for(int i = 0; i < N; i++){
                printf("derivadas[%d] = %lf\n", i, x[i]);
            }
            printf("\n");
            printf("x = %lf\n", punto);
            printf("x0 = %lf\n", x0);


            P = interpolacion_taylor(N, x, punto, x0);

            printf("\n\nEl valor de la interpolacion es: %.12lf\n", P);

            break;

        case 2:
            // double interpolador_lagrange(int N, double *x, double *y, double punto);
            printf("\n\n--- METODO DE LAGRANGE ---\n\n");

            printf("Ingrese el punto en el que se desea calcular: ");
            scanf("%lf", &punto);
            printf("\n");

            for(int i = 0; i < N; i++){
                printf("\nIngrese el valor de x[%d]: ", i);
                scanf("%lf", &x[i]);
                printf(" \nIngrese el valor de y[%d]: ", i);
                scanf("%lf", &y[i]);
            }

            printf("\n\n");
            // imprimimos los datos
            printf("Los datos ingresados son:\n");
            for(int i = 0; i < N; i++){
                printf("x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
            }
            printf("\nEl punto es %lf\n\n", punto);

            P = interpolador_lagrange(N, x, y, punto);

            printf("El valor de la interpolacion es: %.12lf\n", P);

            break;
        case 3:
            // double interpolacion_neville(int N, double *x, double *y, double **Q, double punto);
            printf("\n\n--- METODO DE NEVILLE ---\n\n");
            printf("Ingrese el punto en el que se desea calcular: ");
            scanf("%lf", &punto);
            printf("\n");

            for(int i = 0; i < N; i++){
                printf("\nIngrese el valor de x[%d]: ", i);
                scanf("%lf", &x[i]);
                printf("\nIngrese el valor de y[%d]: ", i);
                scanf("%lf", &y[i]);
            }
            // imprimimos los datos
            printf("\nLos datos ingresados son:\n");
            for(int i = 0; i < N; i++){
                printf("x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
            }
            printf("\n");
            printf("El punto es %lf\n\n", punto);

            P = interpolacion_neville(N, x, y, tabla, punto);

            printf("El valor de la interpolacion es: %.12lf\n", P);

            // tabla de diferencias divididas
            printf("\n");
            printf("La tabla de diferencias divididas es:\n");
            imprimir_matriz(N, tabla);

            break;
        case 4:
            // double interpolador_newton(int N, double *x, double *y, double punto, double **tabla);
            printf("\n\n--- METODO DE NEWTON ---\n\n");
            printf("Ingrese el punto en el que se desea calcular: ");
            scanf("%lf", &punto);
            printf("\n");

            for(int i = 0; i < N; i++){
                printf("\nIngrese el valor de x[%d]: ", i);
                scanf("%lf", &x[i]);
                printf("\nIngrese el valor de y[%d]: ", i);
                scanf("%lf", &y[i]);
            }

            printf("\n\n");
            // imprimitmos los datos
            printf("Los datos ingresados son:\n");
            for(int i = 0; i < N; i++){
                printf("x[%d] = %lf, y[%d] = %lf\n", i, x[i], i, y[i]);
            }
            printf("\n");
            printf("El punto es %lf\n\n", punto);

            P = interpolador_newton(N, x, y, punto, tabla);

            printf("El valor de la interpolacion es: %.12lf\n", P);

            // tabla de diferencias divididas
            printf("\n");
            printf("La tabla de diferencias divididas es:\n");
            imprimir_matriz(N, tabla);

            break;
        case 5:
            printf("\n\nSaliendo...\n\n");
            break;
        default:
            printf("\n\nOpcion no valida\n\n");
            break;
    }
    
    // liberamos la memoria 
    free(x);
    free(y);
    liberar_matriz(N, tabla);

    return 0;
}