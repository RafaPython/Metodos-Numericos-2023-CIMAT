/*
Tarea 5: Métodos de descomposición LU
Nombre: Rafael Alejandro García Ramírez 
email: rafael.ramirez@cimat.mx
Ejecucion : gcc main.c -o main -lm && ./nombre_del_programa NombreArchivo.txt NumeroDeNodos
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "time.h"
#include "solvers.h"

// argv[1] = nombre del archivo a guardar
// argv[2] = numero de nodos

int main(int argc, char *argv[]) {


    if (argc < 3) {
    printf("Error: Se requieren al menos 2 argumentos en la línea de comandos.\n");
    return 1; 
    }

    // nombre del archivo
    char nombre_archivo[strlen(argv[1]) + 1]; // +1 para el carácter nulo '\0'
    strcpy(nombre_archivo, argv[1]);


    // constantes
    double Q = 1000.0, k = 1.0, phi0 = 0.0, phiN = 100.0, l = 10.0;
    int nodos = strtod(argv[2], NULL);
    
    
    nodos -= 1;

    // Abrimos espacio en memoria
    double ** matriz = generar_matriz(nodos);
    double **L = generar_matriz(nodos);
    double **D = generar_matriz(nodos);
    double **L_t = generar_matriz(nodos);

    double *PHI = generar_vector(nodos);
    double *x = generar_vector(nodos);

    double *w = generar_vector(nodos);
    double *v = generar_vector(nodos);

    construir_matriz_calor(nodos, matriz);
    construir_vector_phi(nodos, Q, k, phi0, phiN, l, PHI);

    // cholesky 
    cholesky_LDL(nodos , matriz, L, D);

    // trasnponemos L
    trasponer_matriz(nodos, L, L_t);

    // resolver Lv = PHI
    Lx_b(nodos, L, PHI, v);
    // resolver Dw = v
    Dx_b(nodos, D, v, w);
    // resolver L_t x = w
    Ux_b(nodos, L_t, w, x);

    // guardar en archivo
    guardar_vector(nodos, x, nombre_archivo);


    // liberamos memoria
    liberar_matriz(nodos, matriz);
    liberar_matriz(nodos, L);
    liberar_matriz(nodos, D);
    liberar_matriz(nodos, L_t);
    
    liberar_vector(PHI);
    liberar_vector(x);
    liberar_vector(w);
    liberar_vector(v);

    
    printf("\n----Datos sistema----\n\n");
    printf("Q = %f\n", Q);
    printf("k = %f\n", k);
    printf("phi0 = %f\n",phi0);
    printf("phiN = %f\n",phiN);
    printf("L = %f\n", l);

    printf("\nProceso ejecutado con exito!! Los datos son:\n\n nodos = %d\n Nombre: %s\n\n",nodos + 1, nombre_archivo);
    
    return 0;
}
