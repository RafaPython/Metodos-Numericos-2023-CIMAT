#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"

int main(){

    printf("\n --- METODO DE HERMITE Y SPLINES--- \n");  

    int opcion;
    printf("\n Elija el metodo: \n");
    printf("1. Hermite \n");
    printf("2. Spline Natural \n");
    printf("3. Spline Fijo (Condicionado) \n");
    printf("4. Salir \n");
    printf("Opcion: ");
    scanf("%d", &opcion);
    
    if(opcion == 4){
        printf("\n Saliendo... \n");
        return 0;
    }
    if(opcion < 1 || opcion > 4){
        printf("\n Opcion invalida \n");
        return 0;
    }


    int N;
    printf("\n Ingrese el numero de puntos: ");
    scanf("%d", &N);


    double *x = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *derivadas = (double *)malloc(N * sizeof(double));

    printf("\n Ingrese los puntos x: \n");
    for(int i = 0; i < N; i++){
        printf("x[%d] = ", i);
        scanf("%lf", &x[i]);
    }

    printf("\n Ingrese los puntos y: \n");
    for(int i = 0; i < N; i++){
        printf("y[%d] = ", i);
        scanf("%lf", &y[i]);
    }
    double punto;
    printf("Ingrese el punto a interpolar: ");
    scanf("%lf", &punto);

    int n_condicionados;
    if(opcion == 3){
        
    }

    double P;
    switch (opcion){
        case 1:
            printf("\n\n --- METODO DE HERMITE --- \n");

            printf("\n Ingrese las derivadas: \n");
            for(int i = 0; i < N; i++){
                printf("derivadas[%d] = ", i);
                scanf("%lf", &derivadas[i]);
            }

            printf("Los datos ingresados son: \n");
            printf("x = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", x[i]);
            }
            printf("]\n");

            printf("y = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", y[i]);
            }
            printf("]\n");

            printf("derivadas = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", derivadas[i]);
            }
            printf("]\n");

            printf("punto = %lf\n", punto);

            P = hermite(N, x, y, derivadas, punto);

            printf("\n\nEl valor de la interpolación de Hermite es: %lf\n", P);
            
            break;
        case 2:
            printf("\n\n --- METODO DE SPLINE NATURAL --- \n");

            printf("Los datos ingresados son: \n");
            printf("x = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", x[i]);
            }
            printf("]\n");

            printf("y = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", y[i]);
            }
            printf("]\n");

            printf("punto = %lf\n", punto);

            P = spline_natural(N, x, y, punto);

            printf("\n\nEl valor de la interpolación de Spline Natural es: %lf\n", P);

            break;
        case 3:
            printf("\n\n --- METODO DE SPLINE FIJO (CONDICIONADO) --- \n");
            printf("Cuantos puntos desea condicionar: ");
            scanf("%d", &n_condicionados);
            double *x_condicionados = (double *)malloc(n_condicionados * sizeof(double));
            double *y_condicionados = (double *)malloc(n_condicionados * sizeof(double));
            printf("Ingrese los puntos x condicionados: \n");
            for(int i = 0; i < n_condicionados; i++){
                printf("x[%d] = ", i);
                scanf("%lf", &x_condicionados[i]);
            }

            printf("\n Ingrese los puntos y: \n");
            for(int i = 0; i < n_condicionados; i++){
                printf("y[%d] = ", i);
                scanf("%lf", &y_condicionados[i]);
            }

            printf("Los datos ingresados son: \n");
            printf("x = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", x[i]);
            }
            printf("]\n");

            printf("y = [");
            for(int i = 0; i < N; i++){
                printf("%lf ", y[i]);
            }
            printf("]\n");

            printf("punto = %lf\n", punto);

            // imprimimos los daotos condicionados
            printf("x_condicionados = [");
            for(int i = 0; i < n_condicionados; i++){
                printf("%lf ", x_condicionados[i]);
            }
            printf("]\n");

            printf("y_condicionados = [");
            for(int i = 0; i < n_condicionados; i++){
                printf("%lf ", y_condicionados[i]);
            }
            printf("]\n");

            P = spline_condicionado(N, x, y, x_condicionados, y_condicionados, punto); 

            printf("\n\nEl valor de la interpolación de Spline Fijo (Condicionado) es: %lf\n", P);

            // liberamos memoria 
            free(x_condicionados);
            free(y_condicionados);

            break;
    }

    free(x);
    free(y);
    free(derivadas);

    return 0; 
}