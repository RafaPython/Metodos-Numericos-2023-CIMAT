/*
Autor : Rafael Alejandro García Ramírez
email: rafael.ramirez@cimat.mx

*/


#ifndef LIB_SOLVERS_H

#define LIB_SOLVERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "time.h"
#include "solvers.h"


//Headears de funciones

void doolittle(int N, double **matriz, double **L, double **U);
void crout(int N, double **matriz, double **L, double **U);
void cholesky(int N, double **matriz, double **L);
void cholesky_LDL(int N, double **matriz, double **L, double **D);
void jacobi(int N, double TOL, double **matriz, double *b, double *x0, int iteracionesMAX);
void gauss_seidel(int N, double TOL, double **matriz, double *b, double *x0, int iteracionesMAX);
double metodo_potencia(int N, double **matriz, double *vector, double tolerancia, int max_iteraciones);
double potencia_inversa(int N, double **matriz, double *vector, double tolerancia, int max_iteraciones);
void jacobi_eigen(int N, double **matriz, double **eigenvectores, double tolerancia, int max_iteraciones);
void potencia_k_valores_mayores(int N, int k, double **matriz, double **vectores, double *eigenvalores, double tolerancia, int iteraciones);
void potencia_k_valores_menores(int N, int k, double **matriz, double **vectores, double *eigenvalores, double tolerancia, int iteraciones);



void construir_matriz_calor(int nodos, double **matriz);
void construir_vector_phi(int nodos, double Q, double k, double phi0, double phiN, double L, double *b);


void Lx_b(int N, double **L, double *b, double *x);
void Ux_b(int N, double **U, double *b,double *x);
void Dx_b(int N, double **D, double *b, double *x);

void imprimir_matriz(int N, double **matriz);
void imprimir_vector(int N, double *vector);

void escritura_matriz(FILE *archivo, int filas, int columnas, double **matriz);
void escritura_vector(FILE *archivo, int filas, double *vector);

void multiplicar_matriz_vector(int N, double **A, double *b, double *c);
void multiplicar_matrices(int N, double **A, double **B, double **C); 

void trasponer_matriz(int N, double **A, double **At);
void inversa_triangular_inferior_remplaza(int N, double **L);
void actualizar_rotacion(int N, double **rotacion, double **eigenvectores, int p, int q) ;
void rotar(int N, double **matriz, double **eigenvectores, double theta, int p, int q);
double producto_punto(int N, double *vector1, double *vector2);

double **generar_matriz(int N);
double *generar_vector(int N);
void liberar_matriz(int N, double **matriz);
void liberar_vector(double *vector);
void guardar_vector(int N, double *vector, char *nombre_archivo);

double supremo_diferencias(int N, double *vector1, double *vector2);
int comprobar_solucion(int N, double **A, double *x, double *b, double TOL);
//Funciones

// Para imprimir matrices
void imprimir_matriz(int N, double **matriz){
    // Imprime una matriz de tamaño N x N
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++) printf("%f ", matriz[i][j]);
		printf("\n");
	}
    printf("\n");
}

// Para imprimir vectores
void imprimir_vector(int N, double *vector){
    // Imprime un vector de tamaño N
    for(int i = 0; i < N; i++) printf("%f ", vector[i]);
    printf("\n");
}

// Para multiplicar matrices

void multiplicar_matrices(int N, double **A, double **B, double **C){
    // Multiplica dos matrices A y B y guarda el resultado en C
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) C[i][j] = 0.0;
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) for(int k = 0; k < N; k++) C[i][j] += A[i][k]*B[k][j];
}

void escritura_matriz(FILE *archivo, int filas, int columnas, double **matriz){
    // Leer los valores de la matriz desde el archivo
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            fscanf(archivo, "%lf", &matriz[i][j]);
        }
    }
    fclose(archivo);
}

void escritura_vector(FILE *archivo, int filas, double *vector){
    // Leer los valores de la matriz desde el archivo
    for (int i = 0; i < filas; i++) {
        fscanf(archivo, "%lf", &vector[i]);
    }
    fclose(archivo);
}
void multiplicar_matriz_vector(int N, double **A, double *b, double *c){
    // Multiplica una matriz A por un vector b y guarda el resultado en c
    for(int i = 0; i < N; i++) c[i] = 0.0;
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) c[i] += A[i][j]*b[j];
}

void trasponer_matriz(int N, double **A, double **At){
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) At[i][j] = A[j][i];
}
double producto_punto(int N, double *vector1, double *vector2){
    double producto = 0.0;
    for(int i = 0; i < N; i++) producto += vector1[i]*vector2[i];
    return producto;
}

// ----------------------------------------------------------------------------------------------------------

void doolittle(int N, double **matriz, double **L, double **U){
    /* Realiza la descomposicion de Doolittle de una matriz A en las matrices triangular inferior L y triangular superior U.
    *
    * @ N       El tamano de la matriz A y las matrices L y U (N x N).
    * @ matriz  La matriz de entrada A que se descompone.
    * @ L       La matriz triangular inferior L resultante.
    * @ U       La matriz triangular superior U resultante.
    *
    * Esta funcion implementa el metodo de descomposición de Doolittle para factorizar una matriz A en L y U
    * de manera que A = L * U. Es importante que la matriz A sea no singular y que su elemento diagonal superior no sea cero.
    * Si A no satisface estas condiciones, la funcion imprimira un mensaje de error con informacion de la iteracion de error
    * y terminara el programa.
    */

    if(fabs(matriz[0][0]) <= 1e-16) {printf("Error: matriz imposible. Error Doolittle 0: A[0][0] = 0\n"); exit(1);}

    for (int d = 0; d < N; d++) for(int d2 = 0; d2 < N; d2++) L[d][d2] = U[d][d2] = 0 ;
    for(int i = 0; i < N; i++) L[i][i] = 1.0;

    U[0][0] = matriz[0][0];

    for (int j = 0; j < N; j++ ){
        U[0][j] = matriz[0][j] / L[0][0];
        L[j][0] = matriz[j][0] / U[0][0];
    }
    double acumulador, acumulador1, acumulador2, acumulador3;
    for (int k = 0; k < N - 1; k++){   // k

        acumulador = 0; 
        for(int x = 0; x < k; x++) acumulador += L[k][x]*U[x][k];
        U[k][k] = matriz[k][k] - acumulador;
        
        if (fabs(U[k][k]*L[k][k]) <= 1e-16) {printf("Error: matriz imposible. Error Doolittle 1: U[%d][%d]*L[%d][%d] = 0\n",k,k,k,k); exit(1); }

        for (int x = k + 1; x < N; x++){  // x
            acumulador1 = acumulador2 = 0 ;

            for (int n = 0; n < k; n++) acumulador1 += L[k][n]*U[n][x];
            for (int n = 0; n < k; n++) acumulador2 += L[x][n]*U[n][k];

            U[k][x] = matriz[k][x] - acumulador1;
            L[x][k] = (matriz[x][k] - acumulador2) / U[k][k];
        }
        acumulador3 = 0; 
        for(int n = 0; n < N - 1; n++) acumulador3 += L[N - 1][n]*U[n][N - 1];
        U[N - 1][N - 1] = matriz[N - 1][N - 1] - acumulador3;

        if (fabs(L[N - 1][N - 1]*U[N - 1][N - 1]) <= 1e-16) {printf("Error: matriz imposible. Error Doolitle 2: L[%d][%d]*U[%d][%d] = 0 \n", N-1, N-1, N-1, N-1); exit(1); }
    }
}
void crout(int N, double **matriz, double **L, double **U) {
    /* Descompone una matriz A en las matrices L y U usando el método de Crout.
     *
     * @ N       El tamaño de la matriz A y las matrices L y U.
     * @ matriz  La matriz de entrada A de tamaño N x N.
     * @ L       La matriz triangular inferior L de tamaño N x N.
     * @ U       La matriz triangular superior U de tamaño N x N.
     *
     * Esta funcion implementa el metodo de descomposición de Crout para factorizar una matriz A en L y U
     * de manera que A = L * U. Es importante que la matriz A sea no singular y que su elemento diagonal superior no sea cero.
     * Si A no satisface estas condiciones, la funcion imprimira un mensaje de error con informacion de la iteracion de error
     * y terminara el programa.
     */

    if(fabs(matriz[0][0]) <= 1e-16) { printf("Error: matriz imposible. Error Crout 0: A[0][0] = 0\n"); exit(1);}

    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) U[i][j] = L[i][j] = 0;
    for (int i = 0; i < N; i++) U[i][i] = 1.0; 

    for (int i = 0; i < N; i++){
        for (int j = i; j < N; j++) {
            double acumulador = 0;
            for (int k = 0; k < i; k++) acumulador += L[j][k] * U[k][i];
            L[j][i] = matriz[j][i] - acumulador;
        }
        if (fabs(L[i][i]) <= 1e-16) {printf("Error: matriz imposible. Error Crout 1 : L[%d][%d] = 0\n", i,i); exit(1); }

        for (int j = i + 1; j < N; j++) {
            double acumulador = 0;
            for (int k = 0; k < i; k++) acumulador += L[i][k] * U[k][j];
            U[i][j] = (matriz[i][j] - acumulador) / L[i][i];
        }
    }
}
void Lx_b(int N, double **L, double *b, double *x) {
    /* Resuelve un sistema de ecuaciones lineales triangular inferior Lx = b para x.
     *
     * @ N       El tamano de la matriz triangular inferior L y los vectores b y x.
     * @ L       La matriz triangular inferior de tamano N x N.
     * @ b       El vector de terminos independientes de tamano N.
     * @ x       El vector de solucion, que se actualizara con la solucion despues de la llamada.
     *
     * Esta funcion utiliza sustitucion hacia adelante para resolver el sistema de ecuaciones lineales.
     * La matriz L debe ser triangular inferior, y el vector x debe contener una estimacion inicial
     * de la solucion, que se actualizara con la solucion final.
     */

    double acumulador;
    for (int i = 0; i < N ; i++) {
        acumulador = 0;
        for (int k = 0; k < i; k++) acumulador += L[i][k] * x[k];
        if (fabs(L[i][i]) <= 1e-16) {printf("Error: matriz imposible para Lx = b --> L[%d]L[%d] = 0\n",i,i); exit(1); }
        x[i] = (b[i] - acumulador) / L[i][i];
    }
}
    void Ux_b(int N, double **U, double *b ,double *x){
        /* Resuelve un sistema de ecuaciones lineales triangular superior Ux = b para x.
        *
        * @ N       El tamano de la matriz triangular superior U y los vectores b y x.
        * @ U       La matriz triangular superior de tamano N x N.
        * @ b       El vector de terminos independientes de tamano N.
        * @ x       El vector de solucion, que se actualizara con la solucion despues de la llamada.
        *
        * Esta funcion utiliza sustitucion hacia atras para resolver el sistema de ecuaciones lineales.
        * La matriz U debe ser triangular superior, y el vector x debe contener una estimacion inicial
        * de la solucion, que se actualizara con la solucion final.
        */

        double acumulador;
        for(int i = N - 1; i >= 0; i--){
            acumulador = 0;
            for (int k =  N - 1; k > i; k--) acumulador += U[i][k]*x[k];
            if(fabs(U[i][i]) <= 1e-16) {printf("Error: matriz imposible para Ux = b -- > U[%d][%d] = 0\n",i,i); exit(1);}
            x[i] = (b[i] - acumulador) / U[i][i];
        }  
    }
void Dx_b(int N, double **D, double *b, double *x){
    // D es una matriz diagonal
    // b es un vector
    // x es un vector
    for (int i = 0; i < N; i++) x[i] = b[i] / D[i][i];
}
void cholesky(int N, double **matriz, double **L) {
    /* Descompone una matriz A en la matriz triangular inferior L usando el método de Cholesky.
     *
     * @ N       El tamaño de la matriz A y la matriz L.
     * @ matriz  La matriz de entrada A de tamaño N x N.
     * @ L       La matriz triangular inferior L de tamaño N x N.
     *
     * Esta función realiza la descomposición de Cholesky de la matriz A en la matriz triangular inferior L.
     * Se revisa si la matriz es simétrica y definida positiva, asi como que su primer valor no sea cero.
     * En caso de cumplir con esto, se imprimirá un mensaje y se terminará el programa.
     */
    if (fabs(matriz[0][0]) <= 1e-16) {
        printf("Error: matriz imposible. Error Cholesky 0: A[0][0] = 0\n");
        exit(1);
    }

    // revisamos si la matriz es simétrica
    double **matriz_t = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) matriz_t[i] = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) matriz_t[i][j] = matriz[j][i];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fabs(matriz[i][j] - matriz_t[i][j]) > 1e-16) {
                printf("Error Cholesky: matriz no simetrica. Primer fallo en --> A[%d][%d] != A[%d][%d]\n", i, j, j, i);
                exit(1);
            }
        }
    }

    // Liberar memoria para matriz_t
    for (int i = 0; i < N; i++) {
        free(matriz_t[i]);
    }
    free(matriz_t);

    // Inicializamos la matriz con ceros
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            L[i][j] = 0.0;
        }
    }

    // comenzamos cholesky
    double suma, dummy;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            suma = 0.0;
            if (i == j) {
                for (int k = 0; k < j; k++) suma += L[j][k] * L[j][k];
                if ((dummy = matriz[j][j] - suma) < 0) {
                    printf("Error: matriz no definida positiva. Error Cholesky 1: A[%d][%d] - suma < 0\n", j, j);
                    exit(1);
                }
                L[j][j] = sqrt(dummy);
            } else {
                for (int k = 0; k < j; k++) suma += L[i][k] * L[j][k];
                if (fabs(L[j][j]) <= 1e-16) {
                    printf("Error: Matriz imposible. Error Cholesky 2: L[%d][%d] = 0\n", j, j);
                    exit(1);
                }
                L[i][j] = (matriz[i][j] - suma) / L[j][j];
            }
        }
    }
}
void cholesky_LDL(int N, double **matriz, double **L, double **D) {
    /* Descompone una matriz A en la matriz triangular inferior L y una matriz diagonal D usando el método de Cholesky.
     *
     * @ N       El tamaño de la matriz A y las matrices L y D.
     * @ matriz  La matriz de entrada A de tamaño N x N.
     * @ L       La matriz triangular inferior L de tamaño N x N.
     * @ D       La matriz diagonal D de tamaño N x N.
     *
     * Esta función realiza la descomposición de Cholesky de la matriz A en las matrices L y D.
     * Se revisa si la matriz es simétrica y definida positiva, así como que su primer valor no sea cero.
     * En caso de cumplir con esto, se imprimirá un mensaje y se terminará el programa.
     */
    if (fabs(matriz[0][0]) <= 1e-16) {
        printf("Error: matriz imposible. Error Cholesky 0: A[0][0] = 0\n");
        exit(1);
    }

    // revisamos si la matriz es simétrica
    double **matriz_t = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) matriz_t[i] = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) matriz_t[i][j] = matriz[j][i];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fabs(matriz[i][j] - matriz_t[i][j]) > 1e-16) {
                printf("Error Cholesky: matriz no simétrica. Primer fallo en --> A[%d][%d] != A[%d][%d]\n", i, j, j, i);
                exit(1);
            }
        }
    }

    // Liberar memoria para matriz_t
    for (int i = 0; i < N; i++) {
        free(matriz_t[i]);
    }
    free(matriz_t);

    // Inicializamos las matrices L y D con ceros
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            L[i][j] = 0.0;
            D[i][j] = 0.0;
        }
    }

    // comenzamos cholesky
    double suma, dummy;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            suma = 0.0;
            if (i == j) {
                for (int k = 0; k < j; k++) suma += L[j][k] * D[k][k] * L[j][k];
                if ((dummy = matriz[j][j] - suma) < 0) {
                    printf("Error: matriz no definida positiva. Error Cholesky 1: A[%d][%d] - suma < 0\n", j, j);
                    exit(1);
                }
                D[j][j] = dummy;
                L[j][j] = 1.0;
            } else {
                for (int k = 0; k < j; k++) suma += L[i][k] * D[k][k] * L[j][k];
                if (fabs(D[j][j]) <= 1e-16) {
                    printf("Error: Matriz imposible. Error Cholesky 2: D[%d][%d] = 0\n", j, j);
                    exit(1);
                }
                L[i][j] = (matriz[i][j] - suma) / (D[j][j] * L[j][j]);
                D[i][j] = 0.0;
            }
        }
    }
}
void jacobi(int N, double TOL, double **matriz, double *b, double *x0, int iteracionesMAX){
    /*  Resuelve el sistema de ecuaciones Ax = b por el metodo de Jacobi
    *   
    *   @ N : tamaño de la matriz
    *   @ TOL : tolerancia
    *   @ matriz : matriz de coeficientes
    *   @ b : vector de terminos independientes
    *   @ x0 : vector inicial
    *   @ iteracionesMAX : numero maximo de iteraciones
    * 
    */

    // revisamsos que la matriz sea diagonal dominante y manejamos el error
    double suma;
    int fin = 0;
    for (int i = 0; i < N; i++){
        if (fin == 1) break;
        suma = 0;
        for (int j = 0; j < N; j++){
            if (i != j) suma += fabs(matriz[i][j]);
        }
        if (suma > fabs(matriz[i][i])){
            int opcion;
            printf("ADVERTENCIA Jacobi: La matriz no es diagonal dominante : A[%d][%d]\n", i,i);
            printf("1) Continuar\n2) Salir\nIntroduzca una opcion: ");
            scanf("%d", &opcion);
            switch (opcion){
                case 1:
                    fin = 1;
                    break;
                case 2:
                    exit(1);
                    break;
                default:
                    printf("Opcion no valida\n");
                    exit(1);
                    break;
            }
        }
    }


    // creamos el vector x1
    double *x1 = (double *)malloc(N*sizeof(double));
    for (int i = 0; i < N; i++) x1[i] = 0.0;

    int dummy_iteraciones = iteracionesMAX;
    // comenzamos con jacobbi
    while (iteracionesMAX--){  

        for (int i = 0; i < N; i++){
            suma = 0;
            for (int j = 0; j < N; j++) if(i != j) suma += matriz[i][j]*x0[j];
            x1[i] = (b[i] - suma)/matriz[i][i];
        }

        // revisamos si ya converge
        if (supremo_diferencias(N, x0, x1) < TOL){
            printf("\n\nJacobi: Converge en %d iteraciones\n",dummy_iteraciones - iteracionesMAX -1);
            liberar_vector(x1);
            return;
        }
        // reescribimos x1
        for (int i = 0; i < N; i++) x0[i] = x1[i];
    }
    printf("Error Jacobi 1 : Maximo iteraciones alcanzada : %d\n", dummy_iteraciones);
    // liberamos memoria
    liberar_vector(x1);
}
void gauss_seidel(int N, double TOL, double **matriz, double *b, double *x0, int iteracionesMAX){
    /*  Resuelve el sistema de ecuaciones Ax = b por el metodo de Gauss-Seidel.
    *   
    *   @ N : tamaño de la matriz
    *   @ TOL : tolerancia
    *   @ matriz : matriz de coeficientes
    *   @ b : vector de terminos independientes
    *   @ x0 : vector inicial
    *   @ iteracionesMAX : numero maximo de iteraciones
    * 
    */
    // revisamsos que la matriz sea diagonal dominante y manejamos el error
    double suma;
    int fin = 0;
    for (int i = 0; i < N; i++){
        if (fin == 1) break;
        suma = 0;
        for (int j = 0; j < N; j++){
            if (i != j) suma += fabs(matriz[i][j]);
        }
        if (suma > fabs(matriz[i][i])){
            int opcion;
            printf("ADVERTENCIA Jacobi: La matriz no es diagonal dominante : A[%d][%d]\n", i,i);
            printf("1) Continuar\n2) Salir\nIntroduzca una opcion: ");
            scanf("%d", &opcion);
            switch (opcion){
                case 1:
                    fin = 1;
                    break;
                case 2:
                    exit(1);
                    break;
                default:
                    printf("Opcion no valida\n");
                    exit(1);
                    break;
            }
        }
    }

    double *x1 = (double *) malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) x1[i] = 0.0;

    int dummy_iteraciones = iteracionesMAX;

    double suma1, suma2;
    while(iteracionesMAX--){

        for (int i = 0; i < N; i++){
            suma1 = suma2 = 0;
            for (int j = 0; j < i; j++) suma1 += matriz[i][j] * x1[j];
            for (int j = i+1; j < N; j++) suma2 += matriz[i][j] * x0[j];
            x1[i] = (b[i] - suma1 - suma2) / matriz[i][i];
        }

        // revisamos is ya converge
        if (supremo_diferencias(N, x1, x0) < TOL){
            printf("Gauss-Seidel converge en %d iteraciones\n", dummy_iteraciones - iteracionesMAX - 1);
            liberar_vector(x1);
            return;
        }
        // reescribimos x1
        for (int i = 0; i < N; i++) x0[i] = x1[i];

    }
    printf("Gauss-Seidel no converge en %d iteraciones\n", dummy_iteraciones);
    liberar_vector(x1);
}


double metodo_potencia(int N, double **matriz, double *vector, double tolerancia, int max_iteraciones){
    int copia_max_iteraciones = max_iteraciones;
    double *vector_nuevo = (double *)malloc(N * sizeof(double));
    double max;
    while (max_iteraciones--){

        multiplicar_matriz_vector(N, matriz, vector, vector_nuevo);

        // normalizamos el vector 
        max = vector_nuevo[0];
        for ( int i = 1; i < N; i++) if (fabs(vector_nuevo[i]) > max) max = fabs(vector_nuevo[i]);
        for (int i = 0; i < N; i++) vector_nuevo[i] /= max;

        if (supremo_diferencias(N, vector_nuevo, vector) < tolerancia) {
            free(vector_nuevo);
            printf("Metodo Potencia: Se ha encontrado el eigenvalor más grande después de %d iteraciones.\n", copia_max_iteraciones - max_iteraciones - 1);
            return max;
        }
        // reescribimos el vector
        for (int i = 0; i < N; i++) vector[i] = vector_nuevo[i];
    }
    free(vector_nuevo);
    printf("Error metodo Potencia: No se ha encontrado el eigenvalor más grande después de %d iteraciones.\n", copia_max_iteraciones);
    return -2147483647.0;
}



void inversa_triangular_inferior_remplaza(int N, double **L){
    // revisamos que la matriz sea triangular inferior
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++)
     if (L[i][j] != 0) {
        printf("Error matriz L no triangular inferior : L[%d][%d]\n",i,j);
        return;
    } else if (L[i][i] == 0) {
        printf("Error matriz L no invertible : L[%d][%d]\n",i,i);
        return;
    }
    double suma;
    for(int i = 0; i < N; i++){
        L[i][i] = 1 / L[i][i];
        for (int j = i + 1; j < N; j++){
            suma = 0;
            for (int k = 0; k < j; k++) suma += L[j][k] * L[k][i];
            L[j][i] = -suma / L[j][j];
        }
    }
}



double potencia_inversa(int N, double **matriz, double *vector, double tolerancia, int max_iteraciones){
    int copia_max_iteraciones = max_iteraciones;
    
    double **U = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) U[i] = (double *)malloc(N * sizeof(double));

    double **L = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) L[i] = (double *)malloc(N * sizeof(double));

    // descomponemos la matriz en LU mediante Doolittle
    doolittle(N, matriz, L, U);

    // reescribimos la matriz original por su inversa
    inversa_triangular_inferior_remplaza(N, L);

    // procedemos a realizar el calculo iterativo sobre  Ux^(k+1) = Lx^(k)
    double *dummy = (double *)malloc(N * sizeof(double));
    double *vector_nuevo = (double *)malloc(N * sizeof(double));
    double max;

    while(max_iteraciones--){
        // primero calculamos Lx^(k)
        multiplicar_matriz_vector(N, L, vector, dummy);
        // ahora calculamos Ux^(k+1) = Lx^(k)
        Ux_b(N, U, dummy, vector_nuevo);

        // normalizamos el vector nuevo 
        max = vector_nuevo[0];
        for ( int i = 1; i < N; i++) if (fabs(vector_nuevo[i]) > max) max = fabs(vector_nuevo[i]);
        for (int i = 0; i < N; i++) vector_nuevo[i] /= max;

        // revisamos si se ha cumplido la tolerancia
        if (supremo_diferencias(N, vector_nuevo, vector) < tolerancia) {
            liberar_matriz(N, U);
            liberar_matriz(N, L);
            free(dummy);
            free(vector_nuevo);
            printf("Metodo Potencia Inversa: Se ha encontrado el eigenvalor más pequeño después de %d iteraciones.\n", copia_max_iteraciones - max_iteraciones - 1);
            return 1 / max;
        }
        // reescribimos el vector
        for (int i = 0; i < N; i++) vector[i] = vector_nuevo[i];
        
    }
    liberar_matriz(N, U);
    liberar_matriz(N, L);
    free(dummy);
    free(vector_nuevo);
    printf("Error metodo Potencia Inversa: No se ha encontrado el eigenvalor más pequeño después de %d iteraciones.\n", copia_max_iteraciones);
    return 2147483647.0;
}

void actualizar_rotacion(int N, double **rotacion, double **eigenvectores, int p, int q) {

    double **eigenvectores_nuevo = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) eigenvectores_nuevo[i] = (double *)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) 
        (i == j) ? (eigenvectores_nuevo[i][j] = 1.0): (eigenvectores_nuevo[i][j] = 0.0);

    eigenvectores_nuevo[p][p] = rotacion[p][p];
    eigenvectores_nuevo[p][q] = rotacion[p][q];
    eigenvectores_nuevo[q][p] = rotacion[q][p];
    eigenvectores_nuevo[q][q] = rotacion[q][q];

    // Realiza la multiplicación en una matriz temporal
    double **resultado_multiplicacion = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) resultado_multiplicacion[i] = (double *)malloc(N * sizeof(double));

    multiplicar_matrices(N, eigenvectores, eigenvectores_nuevo, resultado_multiplicacion);

    // Copia el resultado de la multiplicación a eigenvectores
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            eigenvectores[i][j] = resultado_multiplicacion[i][j];
        }
    }

    // Libera la memoria de la matriz temporal y eigenvectores_nuevo
    liberar_matriz(N, resultado_multiplicacion);
    liberar_matriz(N, eigenvectores_nuevo);
}

void rotar(int N, double **matriz, double **eigenvectores, double theta, int p, int q){ 
    // función para realizar la rotación de una matriz R(theta)^T * A * R(theta)
    // donde R(theta) es la matriz de rotación
    // matriz es la matriz a rotar
    // eigenvectores es la matriz de eigenvectores
    // p y q son los índices de la matriz que queremos rotar
    // N es el tamaño de la matriz
    // calculamos el coseno y el seno del ángulo de rotación

    double coseno = cos(theta);
    double seno = sin(theta);

    // inicializamos la matriz de rotación en la identidad
    double **rotacion = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) rotacion[i] = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) (i == j) ? (rotacion[i][j] = 1.0): (rotacion[i][j] = 0.0);

    // actualizamos la matriz de rotación
    rotacion[p][p] = coseno;
    rotacion[p][q] = -seno;
    rotacion[q][p] = seno;
    rotacion[q][q] = coseno;

    // calculamos la matriz rotada
    double **rotacion_transpuesta = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) rotacion_transpuesta[i] = (double *)malloc(N * sizeof(double));
    trasponer_matriz(N, rotacion, rotacion_transpuesta);

    // creamos dummy para cada multiplicación de matrices
    double **dummy = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) dummy[i] = (double *)malloc(N * sizeof(double));
    
    multiplicar_matrices(N, rotacion_transpuesta, matriz, dummy);

    double **dummy2 = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) dummy2[i] = (double *)malloc(N * sizeof(double));

    multiplicar_matrices(N, dummy, rotacion, dummy2);

    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) matriz[i][j] = dummy2[i][j];

    // actualizamos la matriz de eigenvectores
    actualizar_rotacion(N, rotacion, eigenvectores, p, q);

    // liberamos memoria
    liberar_matriz(N, rotacion);
    liberar_matriz(N, rotacion_transpuesta);
    liberar_matriz(N, dummy);
    liberar_matriz(N, dummy2);

}


void jacobi_eigen(int N, double **matriz, double **eigenvectores, double tolerancia, int max_iteraciones){
    double theta, suma;
    int iteraciones = 0;
    for (int iter = 0; iter < max_iteraciones; iter++){
        // hacemos sweep
        for (int n = 0; n < N; n++){
            for(int m = 0; m < n; m++){
                // calculamos el angulo de rotación
                theta = 0.5 * atan2(2 * matriz[m][n], matriz[m][m] - matriz[n][n]);
                // realizamos la rotación
                rotar(N, matriz, eigenvectores, theta, m, n);
            }
            iteraciones++;
        }
        // revisamos si la suma de los elementos fuera de la diagonal es menor que la tolerancia
        suma = 0.0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) suma += matriz[i][j]*matriz[i][j];
        if (suma < tolerancia) {
            printf("Jacobi: Se han encontrado los eigenvalores y eigenvectores después de %d sweeps y %d iteraciones.\n", iter + 1, iteraciones + 1);
            return;
        }
    }
}

void potencia_k_valores_mayores(int N, int k, double **matriz, double **vectores, double *eigenvalores, double tolerancia, int iteraciones){
    double *dummy1 = (double*)malloc(N*sizeof(double)); 
    double *dummy2 = (double*)malloc(N*sizeof(double)); 
    double *temp = (double*)malloc(N*sizeof(double));   

    for(int i = 0; i <k; i++) eigenvalores[i] = 0.0;

    int eigenvalor_anterior, contador;
    double producto, norma;

    for(int i = 0; i < k; i++){
        if( N > 10){
            printf("Encontrando eigenvalores: %d / %d \n", i, k);
        }

        for(int d = 0; d < N; d++) dummy1[d] = 1.0;

        eigenvalor_anterior =  contador = 0;
        
        while (contador < iteraciones){
            for(int d = 0; d < N; d++) dummy2[d] = 0.0;

            for (int j = 0; j < i; j++){
                producto = producto_punto(N, dummy1, vectores[j]);
                for (int jj = 0; jj < N; jj++ ) dummy1[jj] -= producto*vectores[j][jj];
            }
            multiplicar_matriz_vector(N, matriz, dummy1, temp);  // Almacenar el resultado en temp
            eigenvalores[i] = producto_punto(N, dummy1, temp);

            norma = sqrt(producto_punto(N, temp, temp));
            for (int jj = 0; jj < N; jj++) dummy1[jj] = temp[jj] / norma;

            if (fabs(eigenvalores[i] - eigenvalor_anterior) < tolerancia) break;

            eigenvalor_anterior = eigenvalores[i];
            contador++;
        }

        for(int jj = 0; jj < N; jj++ ) vectores[i][jj] = dummy1[jj];
        
    }

    free(dummy1);
    free(dummy2);
    free(temp);
}


void potencia_k_valores_menores(int N, int k, double **matriz, double **vectores, double *eigenvalores, double tolerancia, int iteraciones) {
    double **U = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) U[i] = (double*)malloc(N * sizeof(double));
    double **L = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) L[i] = (double*)malloc(N * sizeof(double));

    // Descomponemos por Doolittle
    doolittle(N, matriz, L, U);

    double *dummy1 = (double*)malloc(N * sizeof(double));
    double *dummy2 = (double*)malloc(N * sizeof(double));
    double *temp = (double*)malloc(N * sizeof(double));

    for (int i = 0; i < k; i++) eigenvalores[i] = 0.0;

    int eigenvalor_anterior, contador;
    double producto, norma;

    for (int i = 0; i < k; i++) {
        // Inicializa dummy1 con un vector aleatorio
        for (int d = 0; d < N; d++) dummy1[d] = ((double)rand() / RAND_MAX) - 0.5;
        for (int d = 0; d < N; d++) dummy2[d] = 0.0;

        eigenvalor_anterior = contador = 0;

        while (contador < iteraciones) {
            for (int j = 0; j < i; j++) {
                producto = producto_punto(N, dummy1, vectores[j]);
                for (int jj = 0; jj < N; jj++) dummy1[jj] -= producto * vectores[j][jj];
            }

            // Solucionamos el sistema LUx = dummy1
            Lx_b(N, L, dummy1, temp);
            Ux_b(N, U, temp, dummy2);

            // Calculamos el eigenvalor como el producto punto entre dummy1 y dummy2
            eigenvalores[i] = 1/ producto_punto(N, dummy1, dummy2);

            norma = sqrt(producto_punto(N, dummy2, dummy2));
            for (int jj = 0; jj < N; jj++) dummy1[jj] = dummy2[jj] / norma;

            if (fabs(eigenvalores[i] - eigenvalor_anterior) < tolerancia) break;

            eigenvalor_anterior = eigenvalores[i];
            contador++;
        }

        for (int jj = 0; jj < N; jj++) vectores[i][jj] = dummy1[jj];
    }

    free(dummy1);
    free(dummy2);
    free(temp);
    liberar_matriz(N, U);
    liberar_matriz(N, L);
}







///////////////////////////////////////////////////////////////////////

double supremo_diferencias(int N, double *vector1, double *vector2){
    // calcula la maxima distancia entre dos like-wise elements de dos vectores.
    double max = 0;
    double dummy = 0;
    for (int i = 0; i < N; i++){
        if((dummy = fabs(vector1[i] - vector2[i])) > max) max = dummy;
    }
    return max;
}

void construir_matriz_calor(int nodos, double **matriz) {
    /* Construye una matriz que describe mediante difetencias finitas el sistema de una barra 1D siendo calentada
     *
     * @ nodos   El numero de nodos en la barra.
     * @ matriz  La matriz de entrada A de tamaño N x N.
     *
     * Las condiciones de fontera se asume que son de Dirichlet, es decir, que los extremos de la barra y que 
     * estos serán tomados en cuenta dentro del vector de terminos independientes. La matriz resultante es una 
     * matriz tridiagonal simetrica con 2 en la diagonal y -1 en las diagonales superior e inferior. En todas 
     * las demás partes es cero. 
    
    */
    // Llenamos la matriz con ceros
    for (int i = 0; i < nodos; i++) for (int j = 0; j < nodos - 1; j++) matriz[i][j] = 0.0;
    // Llenamos la matriz con los valores correspondientes
    for (int i = 0; i < nodos ; i++) {
        for (int j = 0; j < nodos; j++){
            if (i == j) matriz[i][j] = 2.0;
            else matriz[i][j] = 0.0;
            if (i == j + 1 || i == j - 1) matriz[i][j] = -1.0;
        }
    }
}

void construir_vector_phi(int nodos, double Q, double k, double phi0, double phiN, double L, double *b){
    // Construye un vector phi para el sistema de calor
    double delta_x = L/(nodos + 1);
    double dato = (Q*delta_x*delta_x) / k;
    for (int i = 0; i < nodos ; i++) b[i] = dato ;
    b[0] += phi0;
    b[nodos - 1] += phiN; 
}

double **generar_matriz(int N){
    // Genera una matriz cuadrada de tamano N
    double **matriz = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) matriz[i] = (double *)malloc(N * sizeof(double));
    return matriz; 
}
double *generar_vector(int N){
    // Genera un vector de tamano N
    double *vector = (double *)malloc(N * sizeof(double));
    return vector;
}

void liberar_matriz(int N, double **matriz){
    // Libera la memoria de una matriz de tamano N
    for (int i = 0; i < N; i++) free(matriz[i]);
    free(matriz);
}
void liberar_vector(double *vector){
    // Libera la memoria de un vector
    free(vector);
}

void guardar_vector(int N, double *vector, char *nombre_archivo){
    // Guarda un vector en un archivo
    FILE *archivo = fopen(nombre_archivo, "w");
    for (int i = 0; i < N; i++) fprintf(archivo, "%f\n", vector[i]);
    fclose(archivo);
}
int comprobar_solucion(int N, double **A, double *x, double* b, double TOL){
    // Comprueba que Ax - b = 0
    int contador = 0;
    double suma;
    for (int i = 0; i < N; i++){
        suma = 0;
        for (int j = 0; j < N; j++) suma += A[i][j]*x[j];
        if (fabs(suma - b[i]) > TOL) contador++;
    }
    return contador;
}
#endif
