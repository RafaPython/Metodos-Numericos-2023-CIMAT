/*
Tarea 4: Métodos de descomposición LU
Nombre: Rafael Alejandro García Ramírez 
email: rafael.ramirez@cimat.mx
Ejecucion : gcc main.c -o main -lm && ./main archivoMatriz.txt archivoVector.txt TOLERANCIA
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "time.h"
#include "solvers.h"

int main(){
    printf("Programa para solución de matrices generadas tridiagonales mediante Cholesky\n");

    int tamano;
    printf("\n\nIntroduzca el tamaño de la matriz tridiagonal que desee generar: ");
    scanf("%d", &tamano);

    double **matriz = generar_matriz(tamano);

    construir_matriz_calor(tamano, matriz);
    double **L = generar_matriz(tamano);

    cholesky(tamano, matriz, L);

    double ** Lt = generar_matriz(tamano);

    // transponemos la matriz L
    for (int i = 0; i < tamano; i++){
        for (int j = 0; j < tamano; j++){
            Lt[i][j] = L[j][i];
        }
    }

    // guardamos las matrices en archivos
    FILE *archivo_L, *archivo_Lt;
    archivo_L = fopen("matriz_cholesky_L.txt", "w");
    archivo_Lt = fopen("matriz_cholesky_Lt.txt", "w");

    if (archivo_L == NULL) {printf("No se pudo abrir el archivo para L"); return 1;}
    if (archivo_Lt == NULL) {printf("No se pude abrir el archivo para L^T"); return 1;}

    for (int i = 0; i < tamano; i++){
        for (int j = 0; j < tamano; j++){
            fprintf(archivo_L, "%lf ", L[i][j]);
            fprintf(archivo_Lt, "%lf ", Lt[i][j]);
        }
        fprintf(archivo_L, "\n");
        fprintf(archivo_Lt, "\n");
    }

    fclose(archivo_L);
    fclose(archivo_Lt);

    printf("\n\nLas matrices L y L^T de Cholesky han sido guardadas dentro de los archivos 'matriz_cholesky_L.txt' y 'matriz_cholesky_Lt.txt' respectivamente.\n");

    return 0;
}
