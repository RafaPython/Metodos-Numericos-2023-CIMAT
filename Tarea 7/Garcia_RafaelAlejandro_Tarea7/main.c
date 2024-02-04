#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "solvers.h"


/*
ARGUMENTOS:
1. Nombre del archivo con la matriz
2. Tolerancia
3. Maximo de iteraciones

*/
int main(int argc, char *argv[]){   

    if (argc != 4){
        printf("ERROR: Debe ingresar 3 argumentos:\n");
        printf("1. Nombre del archivo con la matriz\n");
        printf("2. Tolerancia\n");
        printf("3. Maximo de iteraciones\n");
        return 1;
    }
    // nombre del archivo con la matriz
    char nombre_archivo[strlen(argv[1]) + 1]; // me gusta este formato 
    strcpy(nombre_archivo, argv[1]);

    printf("Archivo matriz --> %s\n", nombre_archivo);
    printf("Tolerancia ------> %s\n", argv[2]);
    printf("Max iteraciones -> %s  (CUIDADO: evitar uso de notacion cientifica)\n", argv[3]);
    

    double tolerancia = atof(argv[2]);  
    int max_iteraciones = atoi(argv[3]);

    

    // guardamos la matriz en un arreglo
    FILE *archivo_matriz;
    int filas_matriz, columnas_matriz;

    archivo_matriz = fopen(nombre_archivo, "r");
    if (archivo_matriz == NULL) { printf("No se pudo abrir el archivo 1.\n"); return 1; }

    fscanf(archivo_matriz, "%d %d", &filas_matriz, &columnas_matriz);

    double **matriz = generar_matriz(filas_matriz);

    printf("Dimensiones -----> [%d,%d]\n\n", filas_matriz, columnas_matriz);

    escritura_matriz(archivo_matriz, filas_matriz, columnas_matriz, matriz);

    

    printf("\n\nMETODO DE POTENCIA, POTENCIA INVERSA, ENCONTRAR N VALORES PROPIOS MAS GRANDES Y METODO DE JACOBI PARA EIGENVALORES Y EIGENVECTORES\n\n");


    if (filas_matriz < 10){

        int opcion; 
        printf("Ingrese el numero de la opcion que desea realizar:\n");
        printf("1. Metodo de la potencia\n");
        printf("2. Metodo de la potencia inversa\n");
        printf("3. Encontrar los N valores propios mas grandes/pequenios\n");
        printf("4. Metodo de Jacobi para encontrar eigenvalores y eigenvectores\n");
        printf("5. Salir\n");
        printf("Ingrese su opcion: ");
        scanf("%d", &opcion);

        double eigenvalor;
        double *vector = (double *)malloc(filas_matriz * sizeof(double));
        for (int i = 0; i < filas_matriz; i++) vector[i] = 1.0;
        double **eigenvectores = generar_matriz(filas_matriz);

        double *eigenvalores = (double *)malloc(filas_matriz * sizeof(double));
        double **vectores = (double **)malloc(filas_matriz * sizeof(double *));
        for (int i = 0; i < filas_matriz; i++) vectores[i] = (double *)malloc(filas_matriz * sizeof(double));


        int opcionN, numero_eigenvalores;

        switch (opcion) {
        case 1:
            printf("\n--- METODO DE LA POTENCIA---\n\n");
            printf("Matriz original:\n\n");
            imprimir_matriz(filas_matriz, matriz);
            printf("\n");
            eigenvalor = metodo_potencia(filas_matriz, matriz, vector, tolerancia, max_iteraciones);
            printf("\nEigenvalor más grande: %lf\n", eigenvalor);
            printf("\nEigenvector asociado:\n\n");
            imprimir_vector(filas_matriz, vector);
            printf("\n");
            
            break;
        case 2: 
            printf("\n--- METODO DE LA POTENCIA INVERSA---\n\n");
            printf("Matriz original:\n\n");
            imprimir_matriz(filas_matriz, matriz);
            printf("\n");
            eigenvalor = potencia_inversa(filas_matriz, matriz, vector, tolerancia, max_iteraciones);
            printf("\nEigenvalor más pequeño: %lf\n", eigenvalor);
            printf("\nEigenvector asociado:\n\n");
            imprimir_vector(filas_matriz, vector);
            printf("\n");
            break;
        case 3:
            printf("\n--- METODO PARA ENCONTRAR LOS N VALORES PROPIOS MAS GRANDES/PEQUENIOS---\n\n");
            printf("Ingrese el numero de la opcion que desea realizar:\n");
            printf("1. Encontrar los N valores propios mas grandes\n");
            printf("2. Encontrar los N valores propios mas pequenios\n");
            printf("Ingrese su opcion: ");
            scanf("%d", &opcionN);

            if(opcionN == 1){
                printf("Opcion 1: Introduzca la cantidad de valores propios mas grandes que desea encontrar: ");
                scanf("%d", &numero_eigenvalores);
                // imprimimos la matriz original
                printf("Matriz original:\n");
                imprimir_matriz(filas_matriz, matriz);
                printf("\n");
                // imprimimos los N valores propios mas grandes
                printf("Los %d valores propios mas grandes son:\n", numero_eigenvalores);
                potencia_k_valores_mayores(filas_matriz, numero_eigenvalores, matriz, vectores, eigenvalores, tolerancia, max_iteraciones);
                for(int i = 0; i < numero_eigenvalores; i++) printf("%lf\n", eigenvalores[i]);
                printf("\n");
            } else {
                printf("Opcion 2: Introduzca la cantidad de valores propios mas pequenios que desea encontrar: ");
                scanf("%d", &numero_eigenvalores);
                // imprimimos la matriz original
                printf("Matriz original:\n");
                imprimir_matriz(filas_matriz, matriz);
                printf("\n");
                // imprimimos los N valores propios mas grandes
                printf("Los %d valores propios mas pequenios son:\n", numero_eigenvalores);
                potencia_k_valores_menores(filas_matriz, numero_eigenvalores, matriz, vectores, eigenvalores, tolerancia, max_iteraciones);
                for(int i = 0; i < numero_eigenvalores; i++) printf("%lf\n", eigenvalores[i]);
                printf("\n");
            }


            break;
        case 4: 
            printf("\n--- METODO DE JACOBI PARA ENCONTRAR EIGENVALORES Y EIGENVECTORES---\n\n");
            printf("Matriz original:\n\n");
            imprimir_matriz(filas_matriz, matriz);
            printf("\n");
            for (int i = 0; i < filas_matriz; i++) for (int j = 0; j < filas_matriz; j++) (i == j) ? (eigenvectores[i][j] = 1.0): (eigenvectores[i][j] = 0.0);
            jacobi_eigen(filas_matriz, matriz, eigenvectores, tolerancia, max_iteraciones);
            printf("\n\nEigenvalores:\n");
            for (int i = 0; i < filas_matriz; i++) printf("%lf\n", matriz[i][i]);
            printf("\n");
            printf("\n\nEigenvectores:\n");
            imprimir_matriz(filas_matriz, eigenvectores);
            printf("\n");
            break;

        case 5: 
            printf("Saliendo...\n");
            break;
        
        default:
            printf("Opcion no valida.\n");
            break;
        }
        free(vector);
        liberar_matriz(filas_matriz, eigenvectores);
        liberar_matriz(filas_matriz, vectores);
        free(eigenvalores);

    } else {

        printf("\n\nLa dimension de la matriz es mayor a 10. Se guardarán la solución de salida en un archivo para los casos de los eigenvectores y los N valores propios y Jacobi.\n\n");

        int opcion; 
        printf("Ingrese el numero de la opcion que desea realizar:\n");
        printf("1. Metodo de la potencia\n");
        printf("2. Metodo de la potencia inversa\n");
        printf("3. Encontrar los N valores propios mas grandes/pequenios\n");
        printf("4. Metodo de Jacobi para encontrar eigenvalores y eigenvectores\n");
        printf("5. Salir\n");
        printf("Ingrese su opcion: ");
        scanf("%d", &opcion);

        double eigenvalor;
        double *vector = (double *)malloc(filas_matriz * sizeof(double));
        for (int i = 0; i < filas_matriz; i++) vector[i] = 1.0;
        double **eigenvectores = generar_matriz(filas_matriz);

        double *eigenvalores = (double *)malloc(filas_matriz * sizeof(double));
        double **vectores = (double **)malloc(filas_matriz * sizeof(double *));
        for (int i = 0; i < filas_matriz; i++) vectores[i] = (double *)malloc(filas_matriz * sizeof(double));


        int opcionN, numero_eigenvalores;


        FILE *archivo_salida;
        char nombre_salida[100];

        
        FILE *archivo_salida_eigenvalores;
        FILE *archivo_salida_eigenvectores;
        char nombre_salida_eigenvalores[100];
        char nombre_salida_eigenvectores[100];

        switch (opcion){
        case 1: 
            printf("Metodo de la potencia\n");
            eigenvalor = metodo_potencia(filas_matriz, matriz, vector, tolerancia, max_iteraciones);
            printf("Eigenvalor más grande: %lf\n", eigenvalor);
            
            
            sprintf(nombre_salida, "eigenvector_potencia_%d_eigenvalor=%.6f.txt", filas_matriz, eigenvalor);

            
            archivo_salida = fopen(nombre_salida, "w");
            if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
            // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
            fprintf(archivo_salida, "%d %d\n", filas_matriz, 1);
            for (int i = 0; i < filas_matriz; i++) fprintf(archivo_salida, "%f\n", vector[i]);
            fclose(archivo_salida);
            printf("\n\nSe guardó la solución en el archivo %s\n", nombre_salida);
            

            break;
        case 2: 
            printf("Metodo de la potencia inversa\n");
            eigenvalor = potencia_inversa(filas_matriz, matriz, vector, tolerancia, max_iteraciones);
            printf("Eigenvalor más pequeño: %lf\n", eigenvalor);

            // char nombre_salida[100];
            sprintf(nombre_salida, "eigenvector_potencia_inversa_%d_eigenvalor=%.6f.txt", filas_matriz, eigenvalor);

            
            archivo_salida = fopen(nombre_salida, "w");
            if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
            // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
            fprintf(archivo_salida, "%d %d\n", filas_matriz, 1);
            for (int i = 0; i < filas_matriz; i++) fprintf(archivo_salida, "%f\n", vector[i]);
            fclose(archivo_salida);
            printf("\n\nSe guardó la solución en el archivo %s\n", nombre_archivo);
            break;
    

        case 3:
            printf("\n --- METODO PARA ENCONTRAR LOS N VALORES PROPIOS MAS GRANDES/PEQUENIOS ---\n\n");
            printf("Ingrese el numero de la opcion que desea realizar:\n");
            printf("1. Encontrar los N valores propios mas grandes\n");
            printf("2. Encontrar los N valores propios mas pequenios\n");
            printf("Ingrese su opcion: ");
            
            scanf("%d", &opcionN);

            if (opcionN == 1){
                printf("Opcion 1: Introduzca la cantidad de valores propios mas grandes que desea encontrar: ");
                scanf("%d", &numero_eigenvalores);
                potencia_k_valores_mayores(filas_matriz, numero_eigenvalores, matriz, vectores, eigenvalores, tolerancia, max_iteraciones);

                
                sprintf(nombre_salida_eigenvalores, "eigenvalores_kMayores_%d.txt", numero_eigenvalores);

                printf("El nombre de salida es %s\n", nombre_salida_eigenvalores);
                archivo_salida = fopen(nombre_salida_eigenvalores, "w");
                if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
                // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
                fprintf(archivo_salida, "%d %d\n", numero_eigenvalores, 1);
                for (int i = 0; i < numero_eigenvalores; i++) fprintf(archivo_salida, "%f\n", eigenvalores[i]);
                
                printf("\n\nSe guardaron los eigenvalores en en el archivo %s\n", nombre_salida_eigenvalores);

                
                sprintf(nombre_salida_eigenvectores, "eigenvectores_kMayores_%d.txt", numero_eigenvalores);

                archivo_salida = fopen(nombre_salida_eigenvectores, "w");
                if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
                
                // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
                fprintf(archivo_salida, "%d %d\n", filas_matriz, numero_eigenvalores);
                for (int i = 0; i < filas_matriz; i++) for (int j = 0; j < numero_eigenvalores; j++) fprintf(archivo_salida, "%f\n", vectores[i][j]);
                
                printf("\n\nSe guardaron los eigenvectores en en el archivo %s\n", nombre_salida_eigenvectores);
            } else {
                printf("Opcion 2: Introduzca la cantidad de valores propios mas pequenios que desea encontrar: ");
                scanf("%d", &numero_eigenvalores);
                potencia_k_valores_menores(filas_matriz, numero_eigenvalores, matriz, vectores, eigenvalores, tolerancia, max_iteraciones);


                sprintf(nombre_salida_eigenvalores, "eigenvalores_kMenores_%d.txt", numero_eigenvalores);

                printf("El nombre de salida es %s\n", nombre_salida_eigenvalores);
                archivo_salida = fopen(nombre_salida_eigenvalores, "w");
                if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
                // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
                fprintf(archivo_salida, "%d %d\n", numero_eigenvalores, 1);
                for (int i = 0; i < numero_eigenvalores; i++) fprintf(archivo_salida, "%f\n", eigenvalores[i]);

                printf("\n\nSe guardaron los eigenvalores en en el archivo %s\n", nombre_salida_eigenvalores);


                sprintf(nombre_salida_eigenvectores, "eigenvectores_kMenores_%d.txt", numero_eigenvalores);
                
                archivo_salida = fopen(nombre_salida_eigenvectores, "w");
                if (archivo_salida == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }

                // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
                fprintf(archivo_salida, "%d %d\n", filas_matriz, numero_eigenvalores);
                for (int i = 0; i < filas_matriz; i++) for (int j = 0; j < numero_eigenvalores; j++) fprintf(archivo_salida, "%f\n", vectores[i][j]);

                printf("\n\nSe guardaron los eigenvectores en en el archivo %s\n", nombre_salida_eigenvectores);

            }

            break;


        case 4:
            printf("Metodo de Jacobi para encontrar eigenvalores y eigenvectores\n");
            for (int i = 0; i < filas_matriz; i++) for (int j = 0; j < filas_matriz; j++) (i == j) ? (eigenvectores[i][j] = 1.0): (eigenvectores[i][j] = 0.0);

            jacobi_eigen(filas_matriz, matriz, eigenvectores, tolerancia, max_iteraciones);
            
            sprintf(nombre_salida_eigenvalores, "eigenvalores_jacobi_%d.txt", filas_matriz);
            
            sprintf(nombre_salida_eigenvectores, "eigenvectores_jacobi_%d.txt", filas_matriz);

            archivo_salida_eigenvalores = fopen(nombre_salida_eigenvalores, "w");
            if (archivo_salida_eigenvalores == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
            // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
            fprintf(archivo_salida_eigenvalores, "%d %d\n", filas_matriz, 1);
            for (int i = 0; i < filas_matriz; i++) fprintf(archivo_salida_eigenvalores, "%f\n", matriz[i][i]);
            fclose(archivo_salida_eigenvalores);
            printf("\n\nSe guardaron los eigenvalores en en el archivo %s\n", nombre_salida_eigenvalores);
            

            
            archivo_salida_eigenvectores = fopen(nombre_salida_eigenvectores, "w");
            if (archivo_salida_eigenvectores == NULL) { printf("No se pudo abrir el archivo de salida.\n"); return 1; }
            // imprimimos dentro del archivo primero las dimensiones del vector y luego el vector
            fprintf(archivo_salida_eigenvectores, "%d %d\n", filas_matriz, filas_matriz);
            for (int i = 0; i < filas_matriz; i++) for (int j = 0; j < filas_matriz; j++) fprintf(archivo_salida_eigenvectores, "%f\n", eigenvectores[i][j]);
            fclose(archivo_salida_eigenvectores);
            printf("\n\nSe guardaron los eigenvectores en en el archivo %s\n", nombre_salida_eigenvectores);
            break;

        case 5: 
            printf("Saliendo...\n");
            break;

        default:
            printf("Opcion no valida.\n");
            break;
        }
        free(vector);
        liberar_matriz(filas_matriz, eigenvectores);
        liberar_matriz(filas_matriz, vectores);
        free(eigenvalores);
        
    }

    liberar_matriz(filas_matriz, matriz);
    

    return 0;
}