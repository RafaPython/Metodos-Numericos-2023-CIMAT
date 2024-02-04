#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"

/*
1. Nombre archivo matriz 
2. Tolerancia
3. Maximo de iteraciones


NOTA: este programa es solamente para entrega, al momento del examen se usarán archivos para
cada método individualmente. Esto debido a que se maneja memoria para cada caso.

*/

int main(int argc, char *argv[]){

    if (argc != 4){
        printf("Error: Numero de argumentos invalido\n");
        return 1;
    }

    // Nombre del archivo
    char nombre_archivo_matriz[strlen(argv[1]) + 1]; // me gusta este formato
    strcpy(nombre_archivo_matriz, argv[1]);

    double tolerancia = strtod(argv[2], NULL);
    int max_iteraciones = strtod(argv[3], NULL);


    
    printf("Archivo matriz --> %s\n", nombre_archivo_matriz);
    printf("Tolerancia ------> %s\n", argv[2]);
    printf("Iteraciones -----> %s\n", argv[3]);

    FILE *archivo_matriz;
    int filas_matriz, columnas_matriz;


    archivo_matriz = fopen(nombre_archivo_matriz, "r");
    if (archivo_matriz == NULL) { printf("No se pudo abrir el archivo 1.\n"); return 1; }

    fscanf(archivo_matriz, "%d %d", &filas_matriz, &columnas_matriz);
    double **matriz = generar_matriz(filas_matriz);

    escritura_matriz(archivo_matriz, filas_matriz, columnas_matriz, matriz);


    int N = filas_matriz;


            // declaramos todos los vectores y matrices que vamos a necesitar

        // Iteracion subespacio
        double **eigenvalores = generar_matriz(N);
        double **eigenvectores = generar_matriz(N);

        // Rayleigh

        double sigma;
        double *auxRayleigh = (double *)malloc(N * sizeof(double));
        int comprobacion = 0;

        // QR
        double **Q = generar_matriz(filas_matriz);
        double **R = generar_matriz(filas_matriz);
        double **aux = generar_matriz(filas_matriz);


        // Gradiente Conjugado 
        FILE *archivo_vector;
        int filas_vector, columnas_vector;
        char nombre_archivo_vector[100];
        double *b = (double *)malloc(N * sizeof(double));
        double *vector = (double *)malloc(N * sizeof(double));

        // Gradiente Conjugado precondicionado

        double **precondicionador;


    printf("\n\n----RESOLUCIÓN ITERACION SUBESPACIO, RALEIGH, FACTORIZACION QR, GRADIENDE CONJUGADO Y GRADIENTE CONJUGADO PRECONDICIONADO----\n\n");

    if (N < 10){
        int opcion;
        printf("Seleccione el método de resolución:\n");
        printf("1. Iteración del subespacio\n");
        printf("2. Metodo de Raleigh\n");
        printf("3. Factorización QR\n");
        printf("4. Gradiente conjugado\n");
        printf("5. Gradiente conjugado precondicionado\n");
        printf("Opción: ");
        scanf("%d", &opcion);


        switch (opcion)
        {
        case 1:
            printf("\n\n----METODO DE ITERACION DEL SUBESPACIO----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}

            for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) eigenvectores[i][j] = (i == j) ? 1.0 : 0.0;

            // imprimimos la matriz
            printf("Matriz:\n");
            imprimir_matriz(N, matriz);
            printf("\n");
            // hacemos el método de iteración del subespacio
            iteracion_subespacio(N, matriz, eigenvalores, eigenvectores, tolerancia, max_iteraciones);

            printf("Eigenvalores:\n");
            imprimir_matriz(N, eigenvalores);
            printf("\n");

            printf("Eigenvectores:\n");
            imprimir_matriz(N, eigenvectores);
            printf("\n");

            break;
        case 2:
            printf("\n\n----METODO DE RALEIGH----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}
            printf("Ingrese el valor de sigma: ");
            scanf("%lf", &sigma);

            for(int i = 0; i < N; i++) vector[i] = 1.0;

            // imprimimos la matriz
            printf("Matriz:\n");
            imprimir_matriz(N, matriz);
            printf("\n");

            // imprimimos la sigma original
            printf("Sigma original: %lf\n", sigma);
            printf("\n");

            // hacer el método de Raleigh
            sigma = rayleigh(N, sigma, matriz, vector, tolerancia, max_iteraciones);

            // imprimimos la solución
            printf("Vector propio:\n");
            imprimir_vector(N, vector);
            printf("\n");

            // imprimimos la sigma final
            printf("Valor propio asociado: %lf\n", sigma);
            printf("\n");

            // verificamos que A * v - sigma * v = 0
            // le restamos A - sigma * I

            for(int i = 0; i < N; i++) matriz[i][i] -= sigma;
            
            // realizamos la multiplicación matriz * vector
            multiplicar_matriz_vector(N, matriz, vector, auxRayleigh);

            for(int i = 0; i < N; i++) if(fabs(auxRayleigh[i]) > tolerancia) { printf("Error de comprobación Ax -sigmax != en: x[%d] = %lf.\n",i, vector[i]); comprobacion++; }
            if(comprobacion == 0) printf("Comprobación exitosa: Ningún elemento (A - sigma)x != 0\n");
            break;

        case 3:
            printf("\n\n----METODO DE FACTORIZACION QR----\n\n");
            
            // Imprimimos la matriz
            printf("Matriz:\n");
            imprimir_matriz(N, matriz);
            printf("\n");

            // realizamos la factorización QR
            QR(N, N, matriz, Q, R);

            // imprimimos las matrices Q y R
            printf("Matriz Q:\n");
            imprimir_matriz(N, Q);
            printf("\n");

            printf("Matriz R:\n");
            imprimir_matriz(N, R);
            printf("\n");

            // Comprobamos que la factorización es correcta
            comprobacion_QR(N, N, matriz, Q, R, tolerancia);
            

            break;

        case 4:
            printf("\n\n----METODO DE GRADIENTE CONJUGADO----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}
            printf("Ingrese el nombre del archivo del vector b: ");
            scanf("%s", nombre_archivo_vector);

            archivo_vector = fopen(nombre_archivo_vector, "r");
            if (archivo_vector == NULL) { printf("No se pudo abrir el archivo 2.\n"); return 1; }

            fscanf(archivo_vector, "%d %d", &filas_vector, &columnas_vector);
            escritura_vector(archivo_vector, filas_vector, b);
            for(int i = 0; i < filas_vector; i++) vector[i] = 1.0;

            // imprimimos la matriz
            printf("Matriz:\n");
            imprimir_matriz(N, matriz);
            printf("\n");

            // imprimimos el vector b
            printf("Vector b:\n");
            imprimir_vector(N, b);
            printf("\n");

            // hacemos el método de gradiente conjugado
            gradiente_conjugado(N, tolerancia, matriz, b, vector);

            // imprimimos la solución
            printf("\nSolución:\n");
            imprimir_vector(N, vector);
            printf("\n");

            // liberamos memoria
            // liberar_matriz(N, matriz);
            // free(b);
            // free(vector);

            break;

        case 5:
            printf("\n\n----METODO DE GRADIENTE CONJUGADO PRECONDICIONADO----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}
            printf("Ingrese el nombre del archivo del vector b: ");
            scanf("%s", nombre_archivo_vector);

            archivo_vector = fopen(nombre_archivo_vector, "r");
            if (archivo_vector == NULL) { printf("No se pudo abrir el archivo 2.\n"); return 1; }

            fscanf(archivo_vector, "%d %d", &filas_vector, &columnas_vector);
            
            if(filas_matriz != filas_vector) { printf("Las dimensiones no son correctas para el vector y matriz.\n"); break;}
            escritura_vector(archivo_vector, filas_vector, b);

            // calculamos el precondicionador utilizando la matriz 
            precondicionador = generar_matriz(filas_matriz);
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) precondicionador[i][j] = (i == j) ? 1.0 / matriz[i][j] : 0.0;

            for(int i = 0; i < filas_vector; i++) vector[i] = 1.0;

            // imprimimos la matriz
            printf("Matriz:\n");
            imprimir_matriz(N, matriz);
            printf("\n");

            // imprimimos el vector b
            printf("Vector b:\n");
            imprimir_vector(N, b);
            printf("\n");

            // imprimimos el precondicionador
            printf("Precondicionador:\n");
            imprimir_matriz(N, precondicionador);
            printf("\n");

            // hacemos el método de gradiente conjugado precondicionado
            gradiente_conjugado_precondicionado(N, tolerancia, matriz, b, vector, precondicionador);

            // imprimimos la solución
            printf("\nSolución:\n");
            imprimir_vector(N, vector);
            printf("\n");

            // liberamos memoria
            // liberar_matriz(N, matriz);
            // free(b);
            // free(vector);
            liberar_matriz(N, precondicionador);

            break;

        default:
            printf("\n\nOpción no válida.\n\n");
            break;
        }
        
            


    } else {

        printf("\n\nLa dimension de la matriz es mayor a 10. Se guardarán la solución de salida en un archivo.\n\n");
        
        int opcion;
        printf("Seleccione el método de resolución:\n");
        printf("1. Iteración del subespacio\n");
        printf("2. Metodo de Raleigh\n");
        printf("3. Factorización QR\n");
        printf("4. Gradiente conjugado\n");
        printf("5. Gradiente conjugado precondicionado\n");
        printf("Opción: ");
        scanf("%d", &opcion);

        // declaramos todos los vectores y matrices que vamos a necesitar
        
        // Iteración Subespcio

        char nombre_archivo_salida_eigenvalores[100];
        char nombre_archivo_salida_eigenvectores[100];


        // Raleigh


        // QR
        char nombre_archivo_salida_Q[100];
        char nombre_archivo_salida_R[100];


        // QR, Gradiente Conjugado y Gradiente Conjugado precondicionado
    
        
        FILE *archivo_salida_vector;
        char nombre_archivo_salida_vector[100];


        switch (opcion)
        {
        case 1:
            printf("\n\n----METODO DE ITERACION DEL SUBESPACIO----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}
            for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) eigenvectores[i][j] = (i == j) ? 1.0 : 0.0;

            // hacemos el método de iteración del subespacio
            iteracion_subespacio(N, matriz, eigenvalores, eigenvectores, tolerancia, max_iteraciones);

            // guardamos los eigenvalores
            printf("Ingrese el nombre del archivo de salida de los eigenvalores: ");
            scanf("%s", nombre_archivo_salida_eigenvalores);
            guardarMatrizEnArchivo(nombre_archivo_salida_eigenvalores, filas_matriz, columnas_matriz, eigenvalores);

            // guardamos los eigenvectores
            printf("Ingrese el nombre del archivo de salida de los eigenvectores: ");
            scanf("%s", nombre_archivo_salida_eigenvectores);
            guardarMatrizEnArchivo(nombre_archivo_salida_eigenvectores, filas_matriz, columnas_matriz, eigenvectores);

            break;
        case 2:
            printf("\n\n----METODO DE RALEIGH----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}
            printf("Ingrese el valor de sigma: ");
            scanf("%lf", &sigma);

            for(int i = 0; i < N; i++) vector[i] = 1.0;

            // hacemos el método de Raleigh
            sigma = rayleigh(N, sigma, matriz, vector, tolerancia, max_iteraciones);

            printf("El valor propio asociado es: %lf\n", sigma);

            // comprobamos que A * v - sigma * v = 0
            // le restamos A - sigma * I
            for(int i = 0; i < N; i++) matriz[i][i] -= sigma;

            // realizamos la multiplicación matriz * vector
            multiplicar_matriz_vector(N, matriz, vector, auxRayleigh);

            for(int i = 0; i < N; i++) if(fabs(auxRayleigh[i]) > tolerancia) { printf("Error de comprobación Ax -sigmax != en: x[%d] = %lf.\n",i, vector[i]); comprobacion++; }
            if(comprobacion == 0) printf("Comprobación exitosa: Ningún elemento (A - sigma)x != 0\n");

            // guardamos el vector propio
            printf("Ingrese el nombre del archivo de salida del vector propio: ");
            scanf("%s", nombre_archivo_salida_vector);

            archivo_salida_vector = fopen(nombre_archivo_salida_vector, "w");
            if (archivo_salida_vector == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }

            fprintf(archivo_salida_vector, "%d 1\n", filas_matriz);
            for(int i = 0; i < filas_matriz; i++) fprintf(archivo_salida_vector, "%lf\n", vector[i]);

            printf("Archivo de salida guardado.\n");

            break;

        case 3:
            printf("\n\n----METODO DE FACTORIZACION QR----\n\n");
            
            QR(N, N, matriz, Q, R);

            // Comprobamos que la factorización es correcta
            comprobacion_QR(N, N, matriz, Q, R, tolerancia);
            

            // guardamos la matriz Q
            printf("Ingrese el nombre del archivo de salida de la matriz Q: ");
            scanf("%s", nombre_archivo_salida_Q);
            guardarMatrizEnArchivo(nombre_archivo_salida_Q, filas_matriz, columnas_matriz, Q);

            printf("Archivo Q guardado.\n");

            // guardamos la matriz R
            printf("Ingrese el nombre del archivo de salida de la matriz R: ");
            scanf("%s", nombre_archivo_salida_R);
            guardarMatrizEnArchivo(nombre_archivo_salida_R, filas_matriz, columnas_matriz, R);

            printf("Archivo R guardado.\n");


            break;

        case 4:
            printf("\n\n----METODO DE GRADIENTE CONJUGADO----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}

            printf("Ingrese el nombre del archivo del vector b: ");
            scanf("%s", nombre_archivo_vector);

            archivo_vector = fopen(nombre_archivo_vector, "r");
            if (archivo_vector == NULL) { printf("No se pudo abrir el archivo 2.\n"); return 1; }

            fscanf(archivo_vector, "%d %d", &filas_vector, &columnas_vector);
            if(filas_matriz != filas_vector) { printf("Las dimensiones no son correctas para el vector y matriz.\n"); break;}

            escritura_vector(archivo_vector, filas_vector, b);

            for(int i = 0; i < filas_vector; i++) vector[i] = 1.0;

            // realizamos el método de gradiente conjugado
            gradiente_conjugado(N, tolerancia, matriz, b, vector);

            printf("Vector solución:\n");

            imprimir_vector(N, vector);
            printf("\n\n\n\n\n");

            // guardamos la solución
            printf("Ingrese el nombre del archivo de salida de la solución: ");
            scanf("%s", nombre_archivo_salida_vector);

            archivo_salida_vector = fopen(nombre_archivo_salida_vector, "w");
            if (archivo_salida_vector == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }

            fprintf(archivo_salida_vector, "%d 1\n", filas_matriz);
            for(int i = 0; i < filas_matriz; i++) fprintf(archivo_salida_vector, "%lf\n", vector[i]);

            printf("Archivo de salida guardado.\n");

            // liberamos memoria

            // liberar_matriz(N, matriz);
            // free(b);
            // free(vector);


            break;

        case 5:
            printf("\n\n----METODO DE GRADIENTE CONJUGADO PRECONDICIONADO----\n\n");
            if (filas_matriz != columnas_matriz) { printf("La matriz no es cuadrada.\n"); break;}

            printf("Ingrese el nombre del archivo del vector b: ");
            scanf("%s", nombre_archivo_vector);

            archivo_vector = fopen(nombre_archivo_vector, "r");
            if (archivo_vector == NULL) { printf("No se pudo abrir el archivo 2.\n"); return 1; }

            fscanf(archivo_vector, "%d %d", &filas_vector, &columnas_vector);
            if(filas_matriz != filas_vector) { printf("Las dimensiones no son correctas para el vector y matriz.\n"); break;}
            escritura_vector(archivo_vector, filas_vector, b);

            // calculamos el precondicionador utilizando la matriz
            precondicionador = generar_matriz(filas_matriz);
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) precondicionador[i][j] = (i == j) ?  1.0 / matriz[i][j] : 0.0;

            for(int i = 0; i < filas_vector; i++) vector[i] = 1.0;

            // realizamos el método de gradiente conjugado precondicionado
            gradiente_conjugado_precondicionado(N, tolerancia, matriz, b, vector, precondicionador);

            // guardamos la solución
            printf("Ingrese el nombre del archivo de salida de la solución: ");
            scanf("%s", nombre_archivo_salida_vector);

            archivo_salida_vector = fopen(nombre_archivo_salida_vector, "w");
            if (archivo_salida_vector == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }

            fprintf(archivo_salida_vector, "%d 1\n", filas_matriz);
            for(int i = 0; i < filas_matriz; i++) fprintf(archivo_salida_vector, "%lf\n", vector[i]);

            printf("Archivo de salida guardado.\n");
            
            // liberamos memoria
            liberar_matriz(N, precondicionador);
            break;

        default:
            printf("\n\nOpción no válida.\n\n");
            break;
        }

    }


    // liberamos memoria
        liberar_matriz(N, matriz);
        // iteración subespacio
        liberar_matriz(N, eigenvalores);
        liberar_matriz(N, eigenvectores);

        // rayleigh

        free(auxRayleigh);

        // QR
        liberar_matriz(N, Q);
        liberar_matriz(N, R);
        liberar_matriz(N, aux);

        // Gradiente Conjugado
            
        free(b);
        free(vector);

        // Gradiente Conjugado Precondicionado
        // liberar_matriz(N, precondicionador);

    








    return 0;


}