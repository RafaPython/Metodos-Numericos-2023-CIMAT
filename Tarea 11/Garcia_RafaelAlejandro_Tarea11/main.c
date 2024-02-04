#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"
double f1(double x){
    return sin(x);
}
double f2(double x){
    return x*x*log(x);
}
double f3(double x){
    return x*x*exp(-x);
}

int main(){
    int opcion;
    printf("\n--- MINIMOS CUADRADOS E INTEGRACION ---\n\n");
    printf("1. Minimos cuadrados\n");
    printf("2. Integracion\n");
    printf("3. Salir\n");
    printf("\nSeleccione una opcion: ");
    scanf("%d", &opcion);

    if (opcion == 1){
        printf("Metodos de minimos cuadrados:\n");
        printf("Por favor elija el tipo de funcion a ajustar:\n");
        printf("1. Polinomio\n");
        printf("2. Trigonometrica\n");
        printf("3. Exponencial\n");
        printf("4. Salir\n");
        printf("\nSeleccione una opcion: ");
        scanf("%d", &opcion);

        if(opcion == 1){

            printf("\nMetodos de minimos cuadrados para polinomios:\n");
            // introducimos el nombre del archivo x y y a leer
            char nombre_archivo[100];
            printf("\nIntroduzca el nombre del archivo con los datos x: ");
            scanf("%s", nombre_archivo);
            // leemos la primera linea del archivo para saber el numero de puntos
            FILE *archivox = fopen(nombre_archivo, "r");
            int N;
            fscanf(archivox, "%d", &N);
            // fclose(archivox);
            // leemos los datos del archivo
            double *x = (double *)malloc(N*sizeof(double));
            double *y = (double *)malloc(N*sizeof(double));

            printf("\nIntroduzca el nombre del archivo con los datos y: ");
            scanf("%s", nombre_archivo);
            FILE *archivoy = fopen(nombre_archivo, "r");
            fscanf(archivoy, "%d", &N);
            // fclose(archivoy);
            // leemos los datos del archivo
        
            escritura_vector(archivox, N,x);
            escritura_vector(archivox, N,y);

            // creamos el vector de coeficientes
            double *coeficientes = (double *)malloc(N*sizeof(double));
            
            // leemos el regulador lambda
            double lambda;
            printf("\nIntroduzca el valor del regulador lambda: ");
            scanf("%lf", &lambda);

            // calculamos los coeficientes
            interpolacion_cuadrados_polinomio(N, x, y, coeficientes, lambda);

            // guardamos los coeficientes en un archivo
            char nombre_archivo_coeficientes[100];
            printf("\nIntroduzca el nombre del archivo donde guardar los coeficientes: ");
            scanf("%s", nombre_archivo_coeficientes);
            guardar_vector(N, coeficientes, nombre_archivo_coeficientes);

            printf("\nCoeficientes guardados en el archivo %s\n", nombre_archivo_coeficientes);

            printf("Deseea interpolar algun valor? (1: Si, 2: No): ");
            int opcion2;
            scanf("%d", &opcion2);
            if(opcion2 == 1){
                printf("\nIntroduzca el valor a interpolar: ");
                double valor;
                scanf("%lf", &valor);
                double resultado = 0;
                for(int i = 0; i < N; i++){
                    resultado += coeficientes[i]*pow(valor, i);
                }
                printf("\nEl valor interpolado es: %lf\n", resultado);
            }
            if(opcion2 == 2){
                printf("\nSaliendo...\n");
                free(x);
                free(y);
                free(coeficientes);
                return 0;
            }

            // liberamos memoria
            free(x);
            free(y);
            free(coeficientes);
            return 0;

        }
        if(opcion == 2){

            printf("\nMetodos de minimos cuadrados para trigonometrica:\n");
            // introducimos el nombre del archivo x y y a leer
            char nombre_archivo[100];
            printf("\nIntroduzca el nombre del archivo con los datos x: ");
            scanf("%s", nombre_archivo);
            // leemos la primera linea del archivo para saber el numero de puntos
            FILE *archivox = fopen(nombre_archivo, "r");
            int N;
            fscanf(archivox, "%d", &N);
            // fclose(archivox);
            // leemos los datos del archivo
            double *x = (double *)malloc(N*sizeof(double));
            double *y = (double *)malloc(N*sizeof(double));

            printf("\nIntroduzca el nombre del archivo con los datos y: ");
            scanf("%s", nombre_archivo);
            FILE *archivoy = fopen(nombre_archivo, "r");
            fscanf(archivoy, "%d", &N);
            // fclose(archivoy);
            // leemos los datos del archivo
        
            escritura_vector(archivox, N,x);
            escritura_vector(archivox, N,y);

            // creamos el vector de coeficientes
            double *coeficientes = (double *)malloc(N*sizeof(double));
            
            // leemos el regulador lambda
            double lambda;
            printf("\nIntroduzca el valor del regulador lambda: ");
            scanf("%lf", &lambda);

            // calculamos los coeficientes
            interpolacion_cuadrados_coseno(N, x, y, coeficientes, lambda);

            // guardamos los coeficientes en un archivo
            char nombre_archivo_coeficientes[100];
            printf("\nIntroduzca el nombre del archivo donde guardar los coeficientes: ");
            scanf("%s", nombre_archivo_coeficientes);
            guardar_vector(N, coeficientes, nombre_archivo_coeficientes);

            printf("\nCoeficientes guardados en el archivo %s\n", nombre_archivo_coeficientes);

            printf("Deseea interpolar algun valor? (1: Si, 2: No): ");
            int opcion2;
            scanf("%d", &opcion2);
            if(opcion2 == 1){
                printf("\nIntroduzca el valor a interpolar: ");
                double valor;
                scanf("%lf", &valor);
                double resultado = 0;
                for(int i = 0; i < N; i++){
                    resultado += coeficientes[i]*pow(valor, i);
                }
                printf("\nEl valor interpolado es: %lf\n", resultado);
            }
            if(opcion2 == 2){
                printf("\nSaliendo...\n");
                free(x);
                free(y);
                free(coeficientes);
                return 0;
            }

            // liberamos memoria
            free(x);
            free(y);
            free(coeficientes);
            return 0;

        }
        if(opcion == 3){

            printf("\nMetodos de minimos cuadrados para exponencial:\n");
            // introducimos el nombre del archivo x y y a leer
            char nombre_archivo[100];
            printf("\nIntroduzca el nombre del archivo con los datos x: ");
            scanf("%s", nombre_archivo);
            // leemos la primera linea del archivo para saber el numero de puntos
            FILE *archivox = fopen(nombre_archivo, "r");
            int N;
            fscanf(archivox, "%d", &N);
            // fclose(archivox);
            // leemos los datos del archivo
            double *x = (double *)malloc(N*sizeof(double));
            double *y = (double *)malloc(N*sizeof(double));

            printf("\nIntroduzca el nombre del archivo con los datos y: ");
            scanf("%s", nombre_archivo);
            FILE *archivoy = fopen(nombre_archivo, "r");
            fscanf(archivoy, "%d", &N);
            // fclose(archivoy);
            // leemos los datos del archivo
        
            escritura_vector(archivox, N,x);
            escritura_vector(archivox, N,y);

            // creamos el vector de coeficientes
            double *coeficientes = (double *)malloc(N*sizeof(double));
            
            // leemos el regulador lambda
            double lambda;
            printf("\nIntroduzca el valor del regulador lambda: ");
            scanf("%lf", &lambda);

            // calculamos los coeficientes
            interpolacion_cuadrados_exponencial(N, x, y, coeficientes, lambda);

            // guardamos los coeficientes en un archivo
            char nombre_archivo_coeficientes[100];
            printf("\nIntroduzca el nombre del archivo donde guardar los coeficientes: ");
            scanf("%s", nombre_archivo_coeficientes);
            guardar_vector(N, coeficientes, nombre_archivo_coeficientes);

            printf("\nCoeficientes guardados en el archivo %s\n", nombre_archivo_coeficientes);

            printf("Deseea interpolar algun valor? (1: Si, 2: No): ");
            int opcion2;
            scanf("%d", &opcion2);
            if(opcion2 == 1){
                printf("\nIntroduzca el valor a interpolar: ");
                double valor;
                scanf("%lf", &valor);
                double resultado = 0;
                for(int i = 0; i < N; i++){
                    resultado += coeficientes[i]*pow(valor, i);
                }
                printf("\nEl valor interpolado es: %lf\n", resultado);
            }
            if(opcion2 == 2){
                printf("\nSaliendo...\n");
                free(x);
                free(y);
                free(coeficientes);
                return 0;
            }

            // liberamos memoria
            free(x);
            free(y);
            free(coeficientes);
            return 0;

        }
        if(opcion == 4){
            printf("\nSaliendo...\n");
            return 0;
        }
        if (opcion != 1 && opcion != 2 && opcion != 3 && opcion != 4){
            printf("\nOpcion no valida\n");
            return 0;
        }

    }
    if (opcion == 2){
        // pedimos que se elija la funcion a integrar
        int opcionI;
        printf("\nMetodos de integracion:\n");
        printf("Por favor elija el tipo de funcion a integrar:\n");
        printf("1. Sin(x)\n");
        printf("2. x^2*log(x)\n");
        printf("3. x^2*exp(-x)\n");
        printf("4. Salir\n");
        printf("\nSeleccione una opcion: ");
        scanf("%d", &opcionI);

        if (opcionI == 1){
            printf("\nMetodos de integracion para sin(x):\n");
            printf("Por favor elija el metodo de integracion:\n");
            printf("1. Newton-Cotes\n");
            printf("2. Cuadratura de Gauss\n");
            printf("3. Salir\n");
            printf("\nSeleccione una opcion: ");
            int opcionII;
            scanf("%d", &opcionII);

            if(opcionII == 1){
                // preguntamos si se desea abierto o cerrado
                int opcionIII;
                printf("\nMetodos de integracion para sin(x):\n");
                printf("Por favor elija el tipo de metodo de integracion:\n");
                printf("1. Cerrado\n");
                printf("2. Abierto\n");
                printf("3. Salir\n");
                printf("\nSeleccione una opcion: ");
                scanf("%d", &opcionIII);

                if(opcionIII == 1){
                    printf("Metodo de integracion cerrado:\n");
                    double a, b;
                    printf("\nIntroduzca el limite inferior de integracion: ");
                    scanf("%lf", &a);
                    printf("\nIntroduzca el limite superior de integracion: ");
                    scanf("%lf", &b);
                // pedimos el numero n para 
                    int n;
                    printf("\nIntroduzca el numero de subintervalos: ");
                    scanf("%d", &n);
                    // calculamos la integral
                    double resultado = newton_cotes_cerrado(f1, a, b, n);
                    printf("\nEl resultado de la integral es: %lf\n", resultado);

                }
                if(opcionIII == 2){
                    printf("Metodo de integracion abierto:\n");
                    double a, b;
                    printf("\nIntroduzca el limite inferior de integracion: ");
                    scanf("%lf", &a);
                    printf("\nIntroduzca el limite superior de integracion: ");
                    scanf("%lf", &b);

                    // pedimos el numero n para
                    int n;
                    printf("\nIntroduzca el numero de subintervalos: ");
                    scanf("%d", &n);
                    // calculamos la integral
                    double resultado = newton_cotes_abierto(f1, a, b, n);
                    printf("\nEl resultado de la integral es: %lf\n", resultado);

                }
                if (opcionIII == 3){
                    printf("\nSaliendo...\n");
                    return 0;
                }
                if (opcionIII != 1 && opcionIII != 2 && opcionIII != 3){
                    printf("\nOpcion no valida\n");
                    return 0;
                }

            }
            if(opcionII == 2){
                printf("Metodo de integracion de cuadratura de Gauss:\n");
                double a, b;
                printf("\nIntroduzca el limite inferior de integracion: ");
                scanf("%lf", &a);
                printf("\nIntroduzca el limite superior de integracion: ");
                scanf("%lf", &b);

                // pedimos el numero n
                int n;
                printf("\nIntroduzca el numero de subintervalos: ");
                scanf("%d", &n);

                // calculamos la integral
                double resultado = cuadratura_gauss(f1, a, b, n);
                printf("\nEl resultado de la integral es: %lf\n", resultado);

            }
            if (opcionII == 3){
                printf("\nSaliendo...\n");
                return 0;
            }
            if (opcionII != 1 && opcionII != 2 && opcionII != 3){
                printf("\nOpcion no valida\n");
                return 0;
            }

        }


        ////////////////////////////////////////////
        if (opcionI == 2){
            printf("\nMetodos de integracion para x^2*log(x):\n");
            printf("Por favor elija el metodo de integracion:\n");
            printf("1. Newton-Cotes\n");
            printf("2. Cuadratura de Gauss\n");
            printf("3. Salir\n");
            printf("\nSeleccione una opcion: ");
            int opcionII;
            scanf("%d", &opcionII);

            if(opcionII == 1){
                // preguntamos si se desea abierto o cerrado
                int opcionIII;
                printf("\nMetodos de integracion para x^2*log(x):\n");
                printf("Por favor elija el tipo de metodo de integracion:\n");
                printf("1. Cerrado\n");
                printf("2. Abierto\n");
                printf("3. Salir\n");
                printf("\nSeleccione una opcion: ");
                scanf("%d", &opcionIII);

                if(opcionIII == 1){
                    printf("Metodo de integracion cerrado:\n");
                    double a, b;
                    printf("\nIntroduzca el limite inferior de integracion: ");
                    scanf("%lf", &a);
                    printf("\nIntroduzca el limite superior de integracion: ");
                    scanf("%lf", &b);
                // pedimos el numero n para 
                    int n;
                    printf("\nIntroduzca el numero de subintervalos: ");
                    scanf("%d", &n);
                    // calculamos la integral
                    double resultado = newton_cotes_cerrado(f2, a, b, n);
                    printf("\nEl resultado de la integral es: %lf\n", resultado);

                }
                if(opcionIII == 2){
                    printf("Metodo de integracion abierto:\n");
                    double a, b;
                    printf("\nIntroduzca el limite inferior de integracion: ");
                    scanf("%lf", &a);
                    printf("\nIntroduzca el limite superior de integracion: ");
                    scanf("%lf", &b);

                    // pedimos el numero n para
                    int n;
                    printf("\nIntroduzca el numero de subintervalos: ");
                    scanf("%d", &n);
                    // calculamos la integral
                    double resultado = newton_cotes_abierto(f2, a, b, n);
                    printf("\nEl resultado de la integral es: %lf\n", resultado);

                }
                if (opcionIII == 3){
                    printf("\nSaliendo...\n");
                    return 0;
                }
                if (opcionIII != 1 && opcionIII != 2 && opcionIII != 3){
                    printf("\nOpcion no valida\n");
                    return 0;
                }

            }
            if(opcionII == 2){
                printf("Metodo de integracion de cuadratura de Gauss:\n");
                double a, b;
                printf("\nIntroduzca el limite inferior de integracion: ");
                scanf("%lf", &a);
                printf("\nIntroduzca el limite superior de integracion: ");
                scanf("%lf", &b);

                // pedimos el numero n
                int n;
                printf("\nIntroduzca el numero de subintervalos: ");
                scanf("%d", &n);

                // calculamos la integral
                double resultado = cuadratura_gauss(f2, a, b, n);
                printf("\nEl resultado de la integral es: %lf\n", resultado);

            }
            if (opcionII == 3){
                printf("\nSaliendo...\n");
                return 0;
            }
            if (opcionII != 1 && opcionII != 2 && opcionII != 3){
                printf("\nOpcion no valida\n");
                return 0;
            }

            


        }
        ////////////////////////
        if (opcionI == 3){
            printf("\nMetodos de integracion para x^2*exp(-x):\n");
            printf("Por favor elija el metodo de integracion:\n");
            printf("1. Newton-Cotes\n");
            printf("2. Cuadratura de Gauss\n");
            printf("3. Salir\n");
            printf("\nSeleccione una opcion: ");
            int opcionII;
            scanf("%d", &opcionII);

            if(opcionII == 1){
                // preguntamos si se desea abierto o cerrado
                int opcionIII;
                printf("\nMetodos de integracion para x^2*exp(-x):\n");
                printf("Por favor elija el tipo de metodo de integracion:\n");
                printf("1. Cerrado\n");
                printf("2. Abierto\n");
                printf("3. Salir\n");
                printf("\nSeleccione una opcion: ");
                scanf("%d", &opcionIII);

                if(opcionIII == 1){
                    printf("Metodo de integracion cerrado:\n");
                    double a, b;
                    printf("\nIntroduzca el limite inferior de integracion: ");
                    scanf("%lf", &a);
                    printf("\nIntroduzca el limite superior de integracion: ");
                    scanf("%lf", &b);
                // pedimos el numero n para 
                    int n;
                    printf("\nIntroduzca el numero de subintervalos: ");
                    scanf("%d", &n);
                    // calculamos la integral
                    double resultado = newton_cotes_cerrado(f3, a, b, n);
                    printf("\nEl resultado de la integral es: %lf\n", resultado);

                }
                if(opcionIII == 2){
                    printf("Metodo de integracion abierto:\n");
                    double a, b;
                    printf("\nIntroduzca el limite inferior de integracion: ");
                    scanf("%lf", &a);
                    printf("\nIntroduzca el limite superior de integracion: ");
                    scanf("%lf", &b);

                    // pedimos el numero n para
                    int n;
                    printf("\nIntroduzca el numero de subintervalos: ");
                    scanf("%d", &n);
                    // calculamos la integral
                    double resultado = newton_cotes_abierto(f3, a, b, n);
                    printf("\nEl resultado de la integral es: %lf\n", resultado);

                }
                if (opcionIII == 3){
                    printf("\nSaliendo...\n");
                    return 0;
                }
                if (opcionIII != 1 && opcionIII != 2 && opcionIII != 3){
                    printf("\nOpcion no valida\n");
                    return 0;
                }

            }
            if(opcionII == 2){
                printf("Metodo de integracion de cuadratura de Gauss:\n");
                double a, b;
                printf("\nIntroduzca el limite inferior de integracion: ");
                scanf("%lf", &a);
                printf("\nIntroduzca el limite superior de integracion: ");
                scanf("%lf", &b);

                // pedimos el numero n
                int n;
                printf("\nIntroduzca el numero de subintervalos: ");
                scanf("%d", &n);

                // calculamos la integral
                double resultado = cuadratura_gauss(f3, a, b, n);
                printf("\nEl resultado de la integral es: %lf\n", resultado);

            }
            if (opcionII == 3){
                printf("\nSaliendo...\n");
                return 0;
            }
            if (opcionII != 1 && opcionII != 2 && opcionII != 3){
                printf("\nOpcion no valida\n");
                return 0;
            }

        }
        if (opcionI == 4){
            printf("\nSaliendo...\n");
            return 0;
        }
        if (opcionI != 1 && opcionI != 2 && opcionI != 3 && opcionI != 4){
            printf("\nOpcion no valida\n");
            return 0;
        }


    }
    if (opcion == 3){
        printf("\nSaliendo...\n");
        return 0;
    }
    if (opcion != 1 && opcion != 2 && opcion != 3){
        printf("\nOpcion no valida\n");
        return 0;
    }



    return 0; 
}