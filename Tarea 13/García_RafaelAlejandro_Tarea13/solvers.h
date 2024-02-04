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
void gradiente_conjugado(int N, double tolerancia, double **matriz, double *b, double *vector);
void gradiente_conjugado_precondicionado(int N, double tolerancia, double **matriz, double *b, double *vector, double **precondicionador);
double rayleigh(int N, double sigma, double **matriz, double *v0, double tolerancia, int iteraciones);
void QR(int N, int M, double **matriz, double **Q, double **R);

double interpolador_newton(int N, double *x, double *y, double punto, double **tabla);
double interpolacion_neville(int N, double *x, double *y, double **Q, double punto);
double interpolador_lagrange(int N, double *x, double *y, double punto);
double interpolacion_taylor(int N, double *derivadas, double x, double x0);

double hermite(int N, double *x, double *y, double *derivadas, double punto);
double spline_natural(int N, double *x, double *y, double punto);
double spline_fijo(int N, double *x, double *y, double *derivadas, double punto);

double cuadratura_gaussiana(double a, double b, double(*func)(double), int n);


double derivada_hacia_adelante(double (*f)(double), double x, double h);
double derivada_hacia_atras(double (*f)(double), double x, double h);
double derivada_centrada(double (*f)(double), double x, double h);
double derivada_tres_puntos_endpoint(double (*f)(double), double x, double h);
double derivada_tres_puntos_midpoint(double (*f)(double), double x, double h);
double derivada_cinco_puntos_midpoint(double (*f)(double), double x, double h);
double derivada_cinco_puntos_endpoint(double (*f)(double), double x, double h);
double segunda_derivada(double (*f)(double), double x, double h);

double derivada_hacia_adelante_datos(double f1, double f2, double h);
double derivada_hacia_atras_datos(double f1, double f2, double h);
double derivada_centrada_datos(double f1, double f2, double h);
double derivada_tres_puntos_endpoint_datos(double f1, double f2, double f3, double h);
double derivada_tres_puntos_midpoint_datos(double f1, double f2, double h);
double derivada_cinco_puntos_midpoint_datos(double f1, double f2, double f3, double f4, double h);
double derivada_cinco_puntos_endpoint_datos(double f1, double f2, double f3, double f4, double f5, double h);
double segunda_derivada_midpoint_datos(double f1, double f2, double f3, double h);

void jacobiano3(double **J, double x, double y, double z, double h, double (*f1)(double,double,double), double (*f2)(double,double,double), double (*f3)(double,double,double));
void jacobiano2(double **J, double x, double y, double h, double (*f1)(double,double), double (*f2)(double,double));

void hessiana3(double **H, double x, double y, double z, double h, double (*f1)(double, double, double), double (*f2)(double, double, double), double (*f3)(double, double, double));
void hessiana2(double **H, double x, double y, double h, double (*f1)(double, double), double (*f2)(double, double));

void punto_fijo_nolineal(int N, double *x0, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double));
void newton_nolineal(int N, double **J, double *F, double *x0, double tolerancia, double h, int iteraciones, double (*f1)(double,double), double (*f2)(double,double));
void broyden(double *x, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double));
void gradiente_conjugado_nolineal(double *x, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double));


double *euler(int N, double a, double b, double y0, double (*f)(double, double));
double *heun(int N, double a, double b, double y0, double (*f)(double, double));
double derivada_parcial(double (*f)(double, double), double x, double y, double h, int variable);
double *taylor_segundo_orden(int N, double a, double b, double y0, double (*f)(double, double));
double *RK4(int N, double a, double b, double y0, double (*f)(double, double));

double *lotka_volterra(double t, double *xy);
double *euler_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *));
double *heun_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *));
double *RK4_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *));
double *derivada_parcial_sistema(double *(*f)(double, double *), double x, double *y, double h, int variable);
double *taylor_segundo_orden_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *));









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
		for(int j = 0; j < N; j++) printf("%.6f ", matriz[i][j]);
        // memset(matriz[i], 0, N*sizeof(double));
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

void guardarMatrizEnArchivo(const char *nombreArchivo, int filas, int columnas, double **matriz) {
    FILE *archivo = fopen(nombreArchivo, "w");
    
    if (archivo == NULL) {
        printf("No se pudo abrir el archivo %s.\n", nombreArchivo);
        return;
    }
    
    // Guarda las dimensiones de la matriz en la primera fila
    fprintf(archivo, "%d %d\n", filas, columnas);
    
    // Guarda los datos de la matriz en el archivo
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            fprintf(archivo, "%f ", matriz[i][j]);
        }
        fprintf(archivo, "\n"); // Nueva línea para cada fila
    }
    
    fclose(archivo);
    printf("Matriz guardada en %s.\n", nombreArchivo);
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

    
    for (int i = 0; i < N ; i++) {
        x[i] = b[i];
        for (int k = 0; k < i; k++) x[i] -= L[i][k] * x[k];
        if (fabs(L[i][i]) <= 1e-16) {printf("Error: matriz imposible para Lx = b --> L[%d]L[%d] = 0\n",i,i); exit(1); }
        x[i] /= L[i][i];
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
void maximo(double **A, int n, int *m1, int *n1){
    double max = -1;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < i; j++){
            if (i != j && fabs(A[i][j]) > max){
                max = fabs(A[i][j]);
                *m1 = i;
                *n1 = j;
            }
        }
    }
}

void jacobi_eigen(int N, double **matriz, double **eigenvectores, double tolerancia, int max_iteraciones){
    double theta, suma;
    // int iteraciones = 0;
    int px, py;
    for (int iter = 0; iter < max_iteraciones; iter++){
        // hacemos sweep
        maximo(matriz, N, &px, &py);
        theta = 0.5 * atan2(2 * matriz[px][py], matriz[px][px] - matriz[py][py]);
        rotar(N, matriz, eigenvectores, theta, px, py);
        // revisamos si la suma de los elementos fuera de la diagonal es menor que la tolerancia
        suma = 0.0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) suma += matriz[i][j]*matriz[i][j];
        if (suma < tolerancia) {
            // printf("Jacobi: Se han encontrado los eigenvalores y eigenvectores después de %d sweeps y %d iteraciones.\n", iter + 1, iteraciones + 1);
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



void gradiente_conjugado(int N, double tolerancia, double **matriz, double *b, double *vector) {
    int iteraciones = 2 * N;
    double *Ax = (double *)malloc(N * sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *r_nuevo = (double *)malloc(N * sizeof(double));
    double *p = (double *)malloc(N * sizeof(double));
    double *paso_Apk = (double *)malloc(N * sizeof(double));

    // Inicializa rk, pk y rk_nuevo a cero
    double ak, bk;

    // Calcula Ax
    multiplicar_matriz_vector(N, matriz, vector, Ax);

    // Calcula r0
    for (int i = 0; i < N; i++) r[i] = b[i] - Ax[i];

    // Calcula la norma de r0
    double max = 0.0;
    for (int i = 0; i < N; i++) if (fabs(r[i]) > max) max = fabs(r[i]);

    // Comprueba la condición de convergencia
    if (max < tolerancia) {
        printf("Gradiente Conjugado: El vector inicial es una solución aproximada\n");
        free(Ax);
        free(r);
        free(r_nuevo);
        free(p);
        free(paso_Apk);
        return;
    }

    // Inicializa p0
    for (int i = 0; i < N; i++) p[i] = r[i];

    for (int iter = 0; iter < iteraciones; iter++) {
        // Calcula ak
        multiplicar_matriz_vector(N, matriz, p, paso_Apk);
        ak = producto_punto(N, r, r) / producto_punto(N, p, paso_Apk);

        // Calcula x_k+1
        for (int i = 0; i < N; i++) vector[i] += ak * p[i];

        // Calcula r_k+1
        for (int i = 0; i < N; i++) r_nuevo[i] = r[i] - ak * paso_Apk[i];

        // Comprueba la condición de convergencia
        if (sqrt(producto_punto(N, r_nuevo, r_nuevo)) < tolerancia) {
            printf("Gradiente conjugado: el método convergió después de %d iteraciones.\n", iter + 1);
            free(Ax);
            free(r);
            free(r_nuevo);
            free(p);
            free(paso_Apk);
            return;
        }

        // Calcula bk
        bk = producto_punto(N, r_nuevo, r_nuevo) / producto_punto(N, r, r);

        // Actualiza los valores para la siguiente iteración
        for (int i = 0; i < N; i++) p[i] = r_nuevo[i] + bk * p[i];

        // Actualiza r
        for (int i = 0; i < N; i++) r[i] = r_nuevo[i];

    }

    printf("Gradiente conjugado: El método no convergió después de %d iteraciones.\n", 2 * iteraciones);

    // Liberar memoria
    free(Ax);
    free(r);
    free(r_nuevo);
    free(p);
    free(paso_Apk);
}


void gradiente_conjugado_precondicionado(int N, double tolerancia, double **matriz, double *b, double *vector, double **precondicionador) {
    int iteraciones = 2 * N;
    double *Ax = (double *)malloc(N * sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *r_nuevo = (double *)malloc(N * sizeof(double));
    double *z = (double *)malloc(N * sizeof(double));
    double *p = (double *)malloc(N * sizeof(double));
    double *paso_Apk = (double *)malloc(N * sizeof(double));

    // Inicializa rk, pk y rk_nuevo a cero
    double ak, bk;

    // Calcula Az
    multiplicar_matriz_vector(N, matriz, vector, Ax);

    // Calcula r0
    for (int i = 0; i < N; i++) r[i] = b[i] - Ax[i];

    // Aplica precondicionamiento: z = M^(-1) * r
    for (int i = 0; i < N; i++) {
        z[i] = r[i] / precondicionador[i][i];
    }

    // Calcula la norma de r0
    double max = 0.0;
    for (int i = 0; i < N; i++) if (fabs(r[i]) > max) max = fabs(r[i]);

    // Comprueba la condición de convergencia
    if (max < tolerancia) {
        printf("Gradiente Conjugado Precondicionado: El vector inicial es una solución aproximada\n");
        free(Ax);
        free(r);
        free(r_nuevo);
        free(z);
        free(p);
        free(paso_Apk);
        return;
    }

    // Inicializa p0
    for (int i = 0; i < N; i++) p[i] = z[i];
    double rz_interno, rz_nuevo_interno;

    for (int iter = 0; iter < iteraciones; iter++) {
        // Calcula ak
        multiplicar_matriz_vector(N, matriz, p, paso_Apk);

        rz_interno = producto_punto(N, z, r);
        ak = rz_interno / producto_punto(N, p, paso_Apk);

        // Calcula x_k+1
        for (int i = 0; i < N; i++) vector[i] += ak * p[i];

        // Calcula r_k+1
        for (int i = 0; i < N; i++) r_nuevo[i] = r[i] - ak * paso_Apk[i];

        // Comprueba la condición de convergencia
        if (sqrt(producto_punto(N, r_nuevo, r_nuevo)) < tolerancia) {
            printf("Gradiente conjugado precondicionado: el método convergió después de %d iteraciones.\n", iter + 1);
            free(Ax);
            free(r);
            free(r_nuevo);
            free(z);
            free(p);
            free(paso_Apk);
            return;
        }

        // Aplica precondicionamiento: z_k+1 = M^(-1) * r_k+1
        for (int i = 0; i < N; i++) {
            z[i] = r_nuevo[i] / precondicionador[i][i];
        }

        // Calcula bk
        rz_nuevo_interno = producto_punto(N, z, r_nuevo);
        bk = rz_nuevo_interno / rz_interno;

        // Actualiza los valores para la siguiente iteración
        for (int i = 0; i < N; i++) p[i] = z[i] + bk * p[i];

        // Actualiza r
        for (int i = 0; i < N; i++) r[i] = r_nuevo[i];
    }

    printf("Gradiente conjugado precondicionado: El método no convergió después de %d iteraciones.\n", 2 * iteraciones);

    // Liberar memoria
    free(Ax);
    free(r);
    free(r_nuevo);
    free(z);
    free(p);
    free(paso_Apk);
}


double rayleigh(int N, double sigma, double **matriz, double *v0, double tolerancia, int iteraciones){

    // abirmos espacio para la matriz B 
    double **B = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) B[i] = (double *)malloc(N * sizeof(double));

    // abrimos espacio para las matrices B = LU 
    double **U = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) U[i] = (double *)malloc(N * sizeof(double));
    
    double **L = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) L[i] = (double *)malloc(N * sizeof(double));

    double *Ax = (double *)malloc(N * sizeof(double *));

    double *v1 = (double *)malloc(N * sizeof(double));
    for(int i = 0; i < N; i++) v1[i] = 1.0;
    double *y = (double *)malloc(N * sizeof(double));
    double *aux = (double *)malloc(N * sizeof(double));

    double norma;

    for (int i = 1; i <= iteraciones; i++){

        // realizamos la resta de la matriz A - sigma * I
        for(int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) 
                B[i][j] = (i == j) ?  matriz[i][j] - sigma : matriz[i][j]; 

        // ahora resolvemos el sistema de ecuaciones 
        doolittle(N, B, L, U);

        Lx_b(N, L, v0, y);
        Ux_b(N, U, y, v1);        

        norma = sqrt(producto_punto(N, v1, v1));

        for(int i = 0; i < N; i++) v1[i] = v1[i] / norma;

        // calculamos abs(v1 - v0)
        
        for(int i = 0; i < N; i++) aux[i] = v1[i] - v0[i];

        norma = sqrt(producto_punto(N, aux, aux));

        if(norma < tolerancia){
            printf("El método de Rayleigh convergió en %d iteraciones.\n", i);
            // printf("El autovalor es: %lf\n", sigma);
            // printf("El autovector es:\n");
            // imprimir_vector(N, v1);
            // copiamos los valores de v1 a v0
            for(int i = 0; i < N; i++) v0[i] = v1[i];
            printf("\n");

            // liberamos la memoria
            liberar_matriz(N, B);
            liberar_matriz(N, L);
            liberar_matriz(N, U);
            free(Ax);
            free(v1);
            free(y);
            free(aux);

            return sigma;
        }

        // actualizamos el vector v0 = v1
        for(int i = 0; i < N; i++) v0[i] = v1[i];

        // actualizamos el sigma
        multiplicar_matriz_vector(N, matriz, v0, Ax);
    
        sigma = producto_punto(N, v0, Ax) ;
    }

    printf("El método de Rayleigh no convergió en %d iteraciones.\n", iteraciones);

    // liberar memoria
    liberar_matriz(N, B);
    liberar_matriz(N, L);
    liberar_matriz(N, U);
    free(Ax);
    free(v1);
    free(y);
    free(aux);
    return 2147483647.0;
}

void QR(int N, int M, double **matriz, double **Q, double **R) {

    // Copiar la matriz original en una nueva matriz
    double **A = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) A[i] = (double *)malloc(M * sizeof(double));
    for(int i = 0; i < N; i++) for(int j = 0; j < M; j++) A[i][j] = matriz[i][j];

    double norma;
    for (int i = 0; i < N; i++) {
        // Calcular la norma de la columna i
        norma = 0.0;
        for (int j = 0; j < M; j++) norma += A[j][i] * A[j][i];
        norma = sqrt(norma);

        // Actualizar la columna i de R y la columna i de Q
        R[i][i] = norma;
        for (int j = 0; j < M; j++) Q[j][i] = A[j][i] / R[i][i];

        // Actualizar el resto de R y Q
        for (int j = i + 1; j < N; j++) {
            R[i][j] = 0.0;
            for (int k = 0; k < M; k++) R[i][j] += Q[k][i] * A[k][j];
            for (int k = 0; k < M; k++) A[k][j] -= R[i][j] * Q[k][i];
        }
    }

    liberar_matriz(N, A);
}




void comprobacion_QR(int N,int M, double **matriz, double **Q, double **R, double tolerancia){
    double **aux1 = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) aux1[i] = (double *)malloc(N* sizeof(double));

    multiplicar_matrices(N, Q, R, aux1);

    int contador = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(fabs(aux1[i][j] - matriz[i][j]) > tolerancia){
                printf("La factorización QR no es correcta en  A[%d][%d]\n",i,j);
                contador++;
            }
        }
    }

    if(contador == 0) printf("La factorización QR es correcta.\n");
    liberar_matriz(N, aux1);

}


void restar_proyeccion(int p, int q, double **eigenvectores, double*vector){
    double dummy;
    for(int i = 0; i < q; i++){
        dummy = 0;
        for(int j = 0; j < p; j++) dummy += eigenvectores[j][i]*vector[j];
        for(int j = 0; j < p; j++) vector[j] -= dummy*eigenvectores[j][i];
    }
}
void normalizar(double *vector, int longitud) {
    double sumaCuadrados = 0;
    
    // Calcula la suma de los cuadrados de los elementos del vector
    for (int i = 0; i < longitud; i++) {
        sumaCuadrados += vector[i] * vector[i];
    }
    
    // Calcula la magnitud
    double magnitud = sqrt(sumaCuadrados);
    
    // Normaliza el vector dividiendo cada elemento por la magnitud
    for (int i = 0; i < longitud; i++) {
        vector[i] /= magnitud;
    }
}

void iteracion_subespacio(int N, double **matriz,double **eigenvalores, double **eigenvectores, double tolerancia, int iteraciones){
    
    double **p = generar_matriz(N);
    double **pt = generar_matriz(N);
    
    double **L = generar_matriz(N);
    double **U = generar_matriz(N);
    double **aux = generar_matriz(N);

    double *v0 = (double *)malloc(N*sizeof(double));
    double *v1 = (double *)malloc(N*sizeof(double));
    double *y = (double *)malloc(N*sizeof(double));

    // factorizamos A = LU
    doolittle(N, matriz, L, U);

    int diagonal, identidad;

    // inicializamos la matriz p
    for(int i = 0; i < N; i++) p[i][i] = 1.0;

    for(int iter = 0; iter < iteraciones; iter++){

        // Realizamos potencia inversa
        for(int i = 0; i < N; i++){ 
            for(int j = 0; j < N; j++) v0[j] = p[j][i];
            Lx_b(N, L, v0, y);
            Ux_b(N, U, y, v1);
            restar_proyeccion(N, i, p, v1);
            normalizar(v1, N);
            for(int j = 0; j < N; j++) p[j][i] = v1[j]; 

        }

        trasponer_matriz(N,p,pt);
        
        multiplicar_matrices(N, matriz, p, aux);
        multiplicar_matrices(N,pt , aux, eigenvalores);

        // verificamos si eigenvalores es diagonal 
        diagonal = 1;
        for (int j = 0; j < N; j++) {
            if (fabs(eigenvalores[j][j]) < tolerancia){
                diagonal = 0;
                break;
            }
            for (int k = 0; k < N; k++) {
                if (j != k){
                    if (fabs(eigenvalores[j][k]) > tolerancia) {
                        diagonal = 0;
                        break;
                    }
                }
                
            }
        }
        if (diagonal == 1) break;

        jacobi_eigen(N,eigenvalores, eigenvectores, tolerancia,iteraciones);
        
        identidad = 1;
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                if (j != k){
                    if (fabs(producto_punto(N, eigenvectores[j], eigenvectores[k])) > tolerancia) {
                        identidad = 0;
                        break;
                    }
                }
            }
        }

        if (identidad == 1) break;
        

        if(iter == iteraciones - 1){
            printf("Error Iteracion subespacio: el método falló después de %d iteraciones",iteraciones);
            liberar_matriz(N, p);
            liberar_matriz(N, pt);
            liberar_matriz(N, L);
            liberar_matriz(N, U);
            liberar_matriz(N, aux);
            free(v0);
            free(v1);
            free(y);

            return;
        }

        multiplicar_matrices(N, matriz, p, aux);
        


        for(int j = 0; j < N; j++) for(int k = 0; k < N; k++) p[j][k] = aux[j][k];
        
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if (i != j){
                if (fabs(eigenvalores[i][j]) <= tolerancia) eigenvalores[i][j] = 0.0;
            }
        }
    }

    liberar_matriz(N, p);
    liberar_matriz(N, pt);
    liberar_matriz(N, L);
    liberar_matriz(N, U);
    liberar_matriz(N, aux);
    free(v0);
    free(v1);
    free(y);
}



///     INTERPOLACION

double interpolacion_taylor(int N, double *derivadas, double x, double x0){
    double diferencia = x - x0;
    double diferencia_previa = 1;
    double factorial = 1;

    double P = derivadas[0];

    for(int i = 1; i < N ; i++){
        diferencia_previa *= diferencia;
        factorial *= i ;
        P += (derivadas[i] * diferencia_previa) / factorial;  
    }

    return P;
}

double interpolador_lagrange(int N, double *x, double *y, double punto){

    double L, P = 0.0;
    for (int i = 0; i < N; i++){
        L = 1;
        for (int j = 0; j < N; j++){
            if (i != j)
                 L *= (punto - x[j])/(x[i] - x[j]);
        }
        P += y[i] * L;
    }
    return P;
}

double interpolacion_neville(int N, double *x, double *y, double **Q, double punto){

    // hacemos la tabla de diferencias divididas cero 
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) Q[i][j] = 0.0;

    // asignamos los valores de la tabla
    for(int i = 0; i < N; i++) Q[i][0] = y[i];

    // calculamos los valores de la tabla
    for(int i = 1; i < N; i++){
        for(int j = 1; j < i + 1; j++){
            Q[i][j] = ((punto - x[i-j])*Q[i][j-1] - (punto - x[i])*Q[i-1][j-1])/(x[i] - x[i-j]);
        }
    }

    return Q[N-1][N-1];
}


double interpolador_newton(int N, double *x, double *y, double punto, double **tabla){
    // hacemos la tabla de diferencias divididas cero 
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) tabla[i][j] = 0.0;


    // inicializamos la tabla de diferencias divididas
    for(int i = 0; i < N; i++) tabla[i][0] = y[i];

    // calculamos las diferencias divididas
    for(int i = 1; i < N; i++){
        for(int j = 1; j < i + 1; j++){
            tabla[i][j] = (tabla[i][j-1] - tabla[i-1][j-1])/(x[i] - x[i-j]);
        }
    }

    // calculamos el polinomio interpolador
    double L, P = 0.0;
    for(int i = 1; i < N; i++){
        L = 1;
        for(int j = 0; j < N; j++){
            if(j <= i - 1) L *= (punto - x[j]);
        }
        P += tabla[i][i] * L;
    }

    return P + tabla[0][0];
}

double hermite(int N, double *x, double *y, double *derivadas, double punto) {
    // abrimos espacio en la memoria 
    double **Q = generar_matriz(N);

    // inicializamos la matriz en ceros
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) Q[i][j] = 0.0;

    // llenamos la matriz con los valores de y y las derivadas
    for(int i = 0; i < N; i++){
        Q[i][0] = y[i];
        Q[i][1] = derivadas[i];
        if(i != 0) Q[i][1] = (Q[i][0] - Q[i - 1][0]) / (x[i] - x[i - 1]);
    }

    // calculamos los coeficientes
    for(int i = 2; i < N; i++)
         for (int j = 2; j <= i; j++)
             Q[i][j] = (Q[i][j - 1] - Q[i - 1][j - 1]) / (x[i] - x[i - j]);
    
    // calculamos el polinomio
    double producto, P = 0.0;
    for(int i = 0; i < N; i++){
        producto = Q[i][i];
        for(int j = 0; j < i; j++) producto *= (punto - x[j]);
        P += producto;
    }

    // liberamos la memoria 
    liberar_matriz(N, Q);
    
    return P;
}
double spline_natural(int N, double *x, double *y, double punto) {
    double *a = (double *)malloc(N * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *c = (double *)malloc(N * sizeof(double));
    double *d = (double *)malloc(N * sizeof(double));

    double *h = (double *)malloc((N - 1) * sizeof(double));
    double *alpha = (double *)malloc((N - 1) * sizeof(double));

    // Calcular diferencias de x y alpha
    for (int i = 0; i < N - 1; i++) {
        h[i] = x[i + 1] - x[i];
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / (i == 0 ? h[i] : h[i - 1])) * (y[i] - y[i - 1]);
    }

    double *l = (double *)malloc(N * sizeof(double));
    double *u = (double *)malloc(N * sizeof(double));
    double *w = (double *)malloc(N * sizeof(double));

    l[0] = 1.0;
    u[0] = 0.0;
    w[0] = 0.0;

    // Calcular coeficientes l, u y w
    for (int i = 1; i < N - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        w[i] = (alpha[i] - h[i - 1] * w[i - 1]) / l[i];
    }

    c[N - 1] = 0.0;
    l[N - 1] = 1.0;

    // Calcular coeficiente c y retroceder para encontrar a, b y d
    for (int i = N - 2; i >= 0; i--) {
        c[i] = w[i] - u[i] * c[i + 1];
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        a[i] = y[i];
    }

    int intervalo = 0;
    for (int i = 0; i < N - 1; i++) {
        if (x[i] <= punto && punto <= x[i + 1]) {
            intervalo = i;
            break;
        }
    }

    double t = punto - x[intervalo];
    double spline = a[intervalo] + b[intervalo] * t + c[intervalo] * t * t + d[intervalo] * t * t * t;

    free(a);
    free(b);
    free(c);
    free(d);
    free(h);
    free(alpha);
    free(l);
    free(u);
    free(w);

    return spline;
}

double spline_condicionado(int N, double *x, double *y, double *puntos_condicionados, double *valores_condicionados, double punto) {
    double *a = (double *)malloc(N * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *c = (double *)malloc(N * sizeof(double));
    double *d = (double *)malloc(N * sizeof(double));

    // Calcular los coeficientes 'h' y 'alpha'
    double *h = (double *)malloc((N - 1) * sizeof(double));
    double *alpha = (double *)malloc((N - 1) * sizeof(double));

    for (int i = 0; i < N - 1; i++) {
        h[i] = x[i + 1] - x[i];
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / (i == 0 ? h[i] : h[i - 1])) * (y[i] - y[i - 1]);
        
        // Aplicar las condiciones en los puntos dados
        for (int j = 0; j < N - 1; j++) {
            if (x[i] == puntos_condicionados[j]) {
                alpha[i] = valores_condicionados[j];
            }
            if (x[i + 1] == puntos_condicionados[j]) {
                alpha[i] = valores_condicionados[j];
            }
        }
    }

    // Resto del código sin cambios
    double *l = (double *)malloc(N * sizeof(double));
    double *u = (double *)malloc(N * sizeof(double));
    double *w = (double *)malloc(N * sizeof(double));

    l[0] = 1.0; 
    u[0] = 0.0;
    w[0] = 0.0; 

    for (int i = 1; i < N - 1; i++) {
        l[i] = 2 * (h[i - 1] + h[i]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        w[i] = (alpha[i] - alpha[i - 1]) - h[i - 1] * w[i - 1];
        w[i] /= l[i];
    }

    l[N - 1] = 1.0; 
    u[N - 1] = 0.0;
    w[N - 1] = 0.0; 

    // Calcular los coeficientes c, b, d y a
    for (int i = N - 2; i >= 0; i--) {
        c[i] = w[i] - u[i] * c[i + 1];
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3.0;
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        a[i] = y[i];
    }

    // Encontrar el intervalo correspondiente al punto
    int intervalo = 0;
    for (int i = 0; i < N - 1; i++) {
        if (x[i] <= punto && punto <= x[i + 1]) {
            intervalo = i;
            break;
        }
    }

    // Calcular el valor del spline en el punto
    double dx = punto - x[intervalo];
    double resultado = a[intervalo] + b[intervalo] * dx + c[intervalo] * dx * dx + d[intervalo] * dx * dx * dx;

    // Liberar la memoria
    free(a);
    free(b);
    free(c);
    free(d);
    free(h);
    free(alpha);
    free(l);
    free(u);
    free(w);

    return resultado;
}



double newton_cotes_cerrado(double(*func)(double), double a, double b, int n){
    if (n < 1 || n > 4) {
       printf("\n---->Error: n debe ser un entero entre 1 y 5\n");
       return 0.0;
    }
    double h = (b-a)/n;
    double resultado;
    switch (n) {
        case 1:
            resultado = 0.5*h*(func(a) + func(b));
            break;
        case 2:
            resultado = (h/3.0)*(func(a) + 4.0*func(a + h) + func(b));
            break;
        case 3:
            resultado = (3.0*h/8.0)*(func(a) + 3.0*func(a + h) + 3.0*func(a + 2.0*h) + func(b));
            break;
        case 4:
            resultado =(2.0*h/45.0)*(7.0*func(a) + 32.0*func(a + h) + 12.0*func(a + 2.0*h) + 32.0*func(a + 3.0*h) + 7.0*func(b));
            break;
    }
    return resultado;
}

double newton_cotes_abierto(double(*func)(double), double a, double b, int n){
    if (n < 0 || n > 3) {
       printf("\n---->Error: n debe ser un entero entre 0 y 3\n");
       return 0.0;
    }
    double h = (b-a)/(n+2.0);
    double resultado;
    switch (n) {
    case 0:
        resultado = 2*h*func(a + h);
        break;
    case 1:
        resultado = (3.0*h/2.0)*(func(a + h) + func(a + 2.0*h));
        break;
    case 2:
        resultado = (4.0*h/3.0)*(2.0*func(a + h) - func(a + 2.0*h) + 2.0*func(a + 3.0*h));
        break;
    case 3:
        resultado = (5.0*h/24.0)*(11.0*func(a + h) + func(a + 2.0*h) + func(a + 3.0*h) + 11.0*func(a + 4.0*h));
        break;
    }
    return resultado;
}
double cuadratura_gaussiana(double a, double b, double(*func)(double), int n){
    if (n < 1 || n > 4) {
        printf("\n---->Error: n debe ser un entero entre 1 y 4\n");
        return 0.0;
    }
    if (n == 1){
        return 2 * func((a + b) / 2)* (b - a) / 2;
    }
    else if (n == 2){
        double coeficientes[] = {1.0, 1.0};
        double raices[] = {0.5773502692, -0.5773502692};
        double resultado = 0.0;
        double cambio_limites;
        for(int i = 0; i < 2; i++) {
            cambio_limites = ((b-a)*raices[i] + (b+a))/2;
            resultado += coeficientes[i]*func(cambio_limites);
        }
        return resultado*(b - a) / 2.0;
    }
    else if (n == 3){
        double coeficientes[] = {0.555555555555556, 0.888888888888889, 0.555555555555556};
        double raices[] = {-0.774596669241483, 0.0, 0.774596669241483};
        double resultado = 0.0;
        double cambio_limites;
        for(int i = 0; i < 3; i++){
            cambio_limites = ((b-a)*raices[i] + (b+a))/2;
            resultado += coeficientes[i]*func(cambio_limites);
        }
        return resultado*(b - a) / 2.0;
    }
    else if (n == 4){
        double coeficientes[] = {0.652145154862546, 0.652145154862546, 0.347854845137454, 0.347854845137454};
        double raices[] = {-0.339981043584856, 0.339981043584856, -0.861136311594053, 0.861136311594053};
        double resultado = 0.0;
        double cambio_limites;
        for(int i = 0; i < 4; i++){
            cambio_limites = ((b-a)*raices[i] + (b+a))/2;
            resultado += coeficientes[i]*func(cambio_limites);
        }
        return resultado*(b - a) / 2.0;
    }
    return 0.0;
}

double sutherland(double m, double t, double x, double s){
    return m*pow(((1000.0*x)/t),1.5)*(t+s)/(1000.0*x+s);
}


void interpolacion_cuadrados_polinomio(int N, double *x, double *y, double *w, double regulador){

    /*
        @ N: numero de puntos
        @ x: vector de abscisas
        @ y: vector de ordenadas
        @ w: vector de pesos (vacío para escribir los pesos calculados)

    */

    // Creamos la matriz de diseño phi  phi_k = sum phi^T * phi 

    double **phi = generar_matriz(N);
    for (int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            phi[i][j] = pow(x[i], j);
    
    // realizamos el producto phi^T * y
    double **phiT = generar_matriz(N);
    trasponer_matriz(N, phi, phiT);

    // realizamos el producto phi^T * phi
    double **phiTphi = generar_matriz(N);
    multiplicar_matrices(N, phiT, phi, phiTphi);
    for(int i = 0; i < N; i++) phiTphi[i][i] += regulador;

    // multiplicamos phi^T * y
    double *phiTy = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, phiTy);
    


    double *d = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, d);

    // Resolvemos el sistema phiTphi * w = d
    // rellenamos w con valores aleatorios
    double **L = generar_matriz(N);
    double **U = generar_matriz(N);
    doolittle(N, phiTphi, L, U);


    double *dummy = (double *)malloc(N*sizeof(double));
    // resolvemos el sistema L * dummy = d
    Lx_b(N, L, d, dummy);
    // resolvemos el sistema U * w = dummy
    Ux_b(N, U, dummy, w);

    // liberamos la memoria
    liberar_matriz(N, phi);
    liberar_matriz(N, phiT);
    liberar_matriz(N, phiTphi);
    free(phiTy);
    free(d);
    free(dummy);
    liberar_matriz(N, L);
    liberar_matriz(N, U);
}

void interpolacion_cuadrados_coseno(int N, double *x, double *y, double *w, double regulador){

    /*
        @ N: numero de puntos
        @ x: vector de abscisas
        @ y: vector de ordenadas
        @ w: vector de pesos (vacío para escribir los pesos calculados)

    */

    // Creamos la matriz de diseño phi  phi_k = sum phi^T * phi 
    double PI = 3.14159265358979323846;
    double **phi = generar_matriz(N);
    for (int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            phi[i][j] = cos(j*x[i]*PI/6.0);
    
    // realizamos el producto phi^T * y
    double **phiT = generar_matriz(N);
    trasponer_matriz(N, phi, phiT);

    // realizamos el producto phi^T * phi
    double **phiTphi = generar_matriz(N);
    multiplicar_matrices(N, phiT, phi, phiTphi);
    for(int i = 0; i < N; i++) phiTphi[i][i] += regulador;

    // multiplicamos phi^T * y
    double *phiTy = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, phiTy);
    


    double *d = (double *)malloc(N*sizeof(double));
    multiplicar_matriz_vector(N, phiT, y, d);

    // Resolvemos el sistema phiTphi * w = d
    // rellenamos w con valores aleatorios
    double **L = generar_matriz(N);
    double **U = generar_matriz(N);
    doolittle(N, phiTphi, L, U);


    double *dummy = (double *)malloc(N*sizeof(double));
    // resolvemos el sistema L * dummy = d
    Lx_b(N, L, d, dummy);
    // resolvemos el sistema U * w = dummy
    Ux_b(N, U, dummy, w);

    // liberamos la memoria
    liberar_matriz(N, phi);
    liberar_matriz(N, phiT);
    liberar_matriz(N, phiTphi);
    free(phiTy);
    free(d);
    free(dummy);
    liberar_matriz(N, L);
    liberar_matriz(N, U);
}

void interpolacion_cuadrados_exponencial(int N, double *x, double *y, double *w, double regulador){
    
        /*
            @ N: numero de puntos
            @ x: vector de abscisas
            @ y: vector de ordenadas
            @ w: vector de pesos (vacío para escribir los pesos calculados)
    
        */

        // Creamos la matriz de diseño phi  phi_k = sum phi^T * phi
        double **phi = generar_matriz(N);
        double r;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++){ 
                r = (x[i] - x[j])*(x[i] - x[j]);
                phi[i][j] = exp(-r);
            }
                
        
        // realizamos el producto phi^T * y
        double **phiT = generar_matriz(N);
        trasponer_matriz(N, phi, phiT);

        // realizamos el producto phi^T * phi
        double **phiTphi = generar_matriz(N);
        multiplicar_matrices(N, phiT, phi, phiTphi);
        for (int i = 0; i < N; i++)
            phiTphi[i][i] += regulador;

        // multiplicamos phi^T * y
        double *phiTy = (double *)malloc(N * sizeof(double));
        multiplicar_matriz_vector(N, phiT, y, phiTy);

        double *d = (double *)malloc(N * sizeof(double));
        multiplicar_matriz_vector(N, phiT, y, d);

        // Resolvemos el sistema phiTphi * w = d
        // rellenamos w con valores aleatorios
        double **L = generar_matriz(N);
        double **U = generar_matriz(N);
        doolittle(N, phiTphi, L, U);

        double *dummy = (double *)malloc(N * sizeof(double));
        // resolvemos el sistema L * dummy = d
        Lx_b(N, L, d, dummy);
        // resolvemos el sistema U * w = dummy
        Ux_b(N, U, dummy, w);

        // liberamos la memoria
        liberar_matriz(N, phi);
        liberar_matriz(N, phiT);
        liberar_matriz(N, phiTphi);
        free(phiTy);
        free(d);
        free(dummy);
        liberar_matriz(N, L);
        liberar_matriz(N, U);
}

/////// ------ DERIVADAS -------- ///////

/// formulas generales 

double derivada_hacia_adelante(double (*f)(double), double x, double h){
    return (f(x+h)-f(x))/h;
}
double derivada_hacia_atras(double (*f)(double), double x, double h){
    return (f(x)-f(x-h))/h;
}
double derivada_centrada(double (*f)(double), double x, double h){
    return (f(x+h)-f(x-h))/(2*h);
}

/// /////////

double derivada_tres_puntos_endpoint(double (*f)(double), double x, double h){
    return (1/(2*h))*(-3*f(x)+4*f(x+h)-f(x+2*h));
}

double derivada_tres_puntos_midpoint(double (*f)(double), double x, double h){
    return (1/(2*h))*(f(x+h)-f(x-h));
}

double derivada_cinco_puntos_midpoint(double (*f)(double), double x, double h){
    return (1/(12*h))*(f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h));
}

double derivada_cinco_puntos_endpoint(double (*f)(double), double x, double h){
    return (1/(12*h))*(-25*f(x) + 48*f(x+h) - 36*f(x+2*h) + 16*f(x+3*h) - 3*f(x+4*h));
}
double segunda_derivada(double (*f)(double), double x, double h){
    return (1/(h*h))*(f(x-h)-2*f(x)+f(x+h));
}
// ahora lo hacemos para datos en lugar de funciones

double derivada_hacia_adelante_datos(double f1, double f2, double h){
    // f1 es el dato en f(x)
    // f2 es el dato en f(x+h)
    return (f2-f1)/h;
}
double derivada_hacia_atras_datos(double f1, double f2, double h){
    // f1 es el dato en f(x-h)
    // f2 es el dato en f(x)
    return (f2-f1)/h;
}

double derivada_centrada_datos(double f1, double f2, double h){
    // f1 es el dato en f(x-h)
    // f2 es el dato en f(x+h)
    return (f2-f1)/(2*h);
}

double derivada_tres_puntos_endpoint_datos(double f1, double f2, double f3, double h){
    // f1 es el dato en f(x)
    // f2 es el dato en f(x+h)
    // f3 es el dato en f(x+2h)
    return (1/(2*h))*(-3*f1+4*f2-f3);
}
double derivada_tres_puntos_midpoint_datos(double f1, double f2, double h){
    // f1 es el dato en f(x+h)
    // f2 es el dato en f(x-h)
    return (1/(2*h))*(f1-f2);
}
double derivada_cinco_puntos_midpoint_datos( double f1, double f2, double f3, double f4, double h){
    // f1 es el dato en f(x-2h)
    // f2 es el dato en f(x-h)
    // f3 es el dato en f(x+h)
    // f4 es el dato en f(x+2h)
    return (1/12*h)*(f1-8*f2+8*f3-f4);
}
double derivada_cinco_puntos_endpoint_datos(double f1, double f2, double f3, double f4, double f5, double h){
    // f1 es el dato en f(x)
    // f2 es el dato en f(x+h)
    // f3 es el dato en f(x+2h)
    // f4 es el dato en f(x+3h)
    // f5 es el dato en f(x+4h)
    return (1/(12*h))*(-25*f1+48*f2-36*f3+16*f4-3*f5);
}
double segunda_derivada_midpoint_datos(double f1, double f2, double f3, double h){
    // f1 es el dato en f(x-h)
    // f2 es el dato en f(x)
    // f3 es el dato en f(x+h)
    return (1/(h*h))*(f1-2*f2+f3);

}


void jacobiano3(double **J, double x, double y, double z, double h, double (*f1)(double,double,double), double (*f2)(double,double,double), double (*f3)(double,double,double)){
    J[0][0] = (1/(12*h))*(-25*f1(x,y,z) + 48*f1(x+h,y,z) - 36*f1(x+2*h,y,z) + 16*f1(x+3*h,y,z) - 3*f1(x+4*h,y,z));
    J[0][1] = (1/(12*h))*(-25*f1(x,y,z) + 48*f1(x,y+h,z) - 36*f1(x,y+2*h,z) + 16*f1(x,y+3*h,z) - 3*f1(x,y+4*h,z));
    J[0][2] = (1/(12*h))*(-25*f1(x,y,z) + 48*f1(x,y,z+h) - 36*f1(x,y,z+2*h) + 16*f1(x,y,z+3*h) - 3*f1(x,y,z+4*h));

    J[1][1] = (1/(12*h))*(-25*f2(x,y,z) + 48*f2(x,y+h,z) - 36*f2(x,y+2*h,z) + 16*f2(x,y+3*h,z) - 3*f2(x,y+4*h,z));
    J[2][2] = (1/(12*h))*(-25*f3(x,y,z) + 48*f3(x,y,z+h) - 36*f3(x,y,z+2*h) + 16*f3(x,y,z+3*h) - 3*f3(x,y,z+4*h));
    J[1][2] = (1/(12*h))*(-25*f2(x,y,z) + 48*f2(x,y,z+h) - 36*f2(x,y,z+2*h) + 16*f2(x,y,z+3*h) - 3*f2(x,y,z+4*h));
    J[2][1] = (1/(12*h))*(-25*f3(x,y,z) + 48*f3(x,y+h,z) - 36*f3(x,y+2*h,z) + 16*f3(x,y+3*h,z) - 3*f3(x,y+4*h,z));
    J[1][0] = (1/(12*h))*(-25*f2(x,y,z) + 48*f2(x+h,y,z) - 36*f2(x+2*h,y,z) + 16*f2(x+3*h,y,z) - 3*f2(x+4*h,y,z));
    J[2][0] = (1/(12*h))*(-25*f3(x,y,z) + 48*f3(x+h,y,z) - 36*f3(x+2*h,y,z) + 16*f3(x+3*h,y,z) - 3*f3(x+4*h,y,z));
}

void jacobiano2(double **J, double x, double y, double h, double (*f1)(double,double), double (*f2)(double,double)){
    J[0][0] = (1/(12*h))*(-25*f1(x,y) + 48*f1(x+h,y) - 36*f1(x+2*h,y) + 16*f1(x+3*h,y) - 3*f1(x+4*h,y));
    J[0][1] = (1/(12*h))*(-25*f1(x,y) + 48*f1(x,y+h) - 36*f1(x,y+2*h) + 16*f1(x,y+3*h) - 3*f1(x,y+4*h));

    J[1][1] = (1/(12*h))*(-25*f2(x,y) + 48*f2(x,y+h) - 36*f2(x,y+2*h) + 16*f2(x,y+3*h) - 3*f2(x,y+4*h));
    J[1][0] = (1/(12*h))*(-25*f2(x,y) + 48*f2(x+h,y) - 36*f2(x+2*h,y) + 16*f2(x+3*h,y) - 3*f2(x+4*h,y));
}

void hessiana3(double **H, double x, double y, double z, double h, double (*f1)(double, double, double), double (*f2)(double, double, double), double (*f3)(double, double, double)){
    H[0][0] = (1 / (h * h)) * (f1(x + h, y, z) - 2 * f1(x, y, z) + f1(x - h, y, z));
    H[0][1] = (1 / (h * h)) * (f1(x + h, y + h, z) - 2 * f1(x, y + h, z) + f1(x - h, y + h, z) - 2 * f1(x + h, y, z) + 4 * f1(x, y, z) - 2 * f1(x - h, y, z) + f1(x + h, y - h, z) - 2 * f1(x, y - h, z) + f1(x - h, y - h, z));
    H[0][2] = (1 / (h * h)) * (f1(x + h, y, z + h) - 2 * f1(x, y, z + h) + f1(x - h, y, z + h) - 2 * f1(x + h, y, z) + 4 * f1(x, y, z) - 2 * f1(x - h, y, z) + f1(x + h, y, z - h) - 2 * f1(x, y, z - h) + f1(x - h, y, z - h));

    H[1][0] = H[0][1];
    H[1][1] = (1 / (h * h)) * (f2(x + h, y, z) - 2 * f2(x, y, z) + f2(x - h, y, z));
    H[1][2] = (1 / (h * h)) * (f2(x, y + h, z) - 2 * f2(x, y, z) + f2(x, y - h, z));

    H[2][0] = H[0][2];
    H[2][1] = H[1][2];
    H[2][2] = (1 / (h * h)) * (f3(x + h, y, z) - 2 * f3(x, y, z) + f3(x - h, y, z));
}
void hessiana2(double **H, double x, double y, double h, double (*f1)(double, double), double (*f2)(double, double)){
    H[0][0] = (1 / (h * h)) * (f1(x + h, y) - 2 * f1(x, y) + f1(x - h, y));
    H[0][1] = (1 / (h * h)) * (f1(x + h, y + h) - 2 * f1(x, y + h) + f1(x - h, y + h) - 2 * f1(x + h, y) + 4 * f1(x, y) - 2 * f1(x - h, y) + f1(x + h, y - h) - 2 * f1(x, y - h) + f1(x - h, y - h));
    
    H[1][0] = H[0][1];
    H[1][1] = (1 / (h * h)) * (f2(x + h, y) - 2 * f2(x, y) + f2(x - h, y));
}

void punto_fijo_nolineal(int N, double *x0, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double)) {
    int i_convergencia = 0;
    double x1_anterior, x2_anterior;
    while (iteraciones--) {
        for (int j = 0; j < N; j++) {
            x1_anterior = x0[0];
            x2_anterior = x0[1];
            x0[0] = f1(x0[0], x0[1]);
            x0[1] = f2(x0[0], x0[1]);
        }
        i_convergencia++;
        if (fabs(x0[0] - x1_anterior) < tolerancia && fabs(x0[1] - x2_anterior) < tolerancia) {
            printf("Punto Fijo no lineal: convergencia alcanzada en %d iteraciones\n", i_convergencia);
            return;
        }
    }
    printf("Punto Fijo no lineal: NO converge en %d iteraciones\n", i_convergencia);
}
void newton_nolineal(int N, double **J, double *F, double *x0, double tolerancia, double h, int iteraciones, double (*f1)(double,double), double (*f2)(double,double)){
    int i_convergencia = 0;
    double *x_anterior = generar_vector(N);
    double dummy;
    while(iteraciones--){
        F[0] = f1(x0[0], x0[1]);
        F[1] = f2(x0[0], x0[1]);
        jacobiano2(J, x0[0], x0[1], h, f1, f2);

        // resolvemos el sistema Jy = -F
        dummy = -J[1][0] / J[0][0];
        J[1][1] -= dummy * J[0][1];
        F[1] -= dummy * F[0];

        x_anterior[1] = F[1] / J[1][1];
        x_anterior[0] = (F[0] - J[0][1] * x_anterior[1]) / J[0][0];

        x0[0] += x_anterior[0];
        x0[1] += x_anterior[1];

        // printf("x1: %lf, x2: %lf\n", x0[0], x0[1]);

        if(fabs(x_anterior[0]) < tolerancia && fabs(x_anterior[1]) < tolerancia){
            printf("Newton no lineal: convergencia alcanzada en %d iteraciones\n", i_convergencia);
            free(x_anterior);
            return;
        }
        i_convergencia++;
    }
    printf("Newton no lineal: no converge en %d iteraciones\n", i_convergencia);
    free(x_anterior);
}

void broyden(double *solucion, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double)) {
    int iteracion = 0;
    double F[2], s[2], u[2], v_anterior[2], w_anterior[2], yk[2], z[2];
    double matriz_jacobiana[2][2];
    double producto_punto;

    matriz_jacobiana[0][0] = 1.0, matriz_jacobiana[0][1] = 1.0;
    matriz_jacobiana[1][0] = 2.0, matriz_jacobiana[1][1] = 2.0;

    while (iteraciones--) {
        w_anterior[0] = v_anterior[0];
        w_anterior[1] = v_anterior[1];
        F[0] = f1(solucion[0], solucion[1]);
        F[1] = f2(solucion[0], solucion[1]);
        v_anterior[0] = F[0];
        v_anterior[1] = F[1];

        yk[0] = v_anterior[0] - w_anterior[0];
        yk[1] = v_anterior[1] - w_anterior[1];

        z[0] = -matriz_jacobiana[0][0] * yk[0] - matriz_jacobiana[0][1] * yk[1];
        z[1] = -matriz_jacobiana[1][0] * yk[0] - matriz_jacobiana[1][1] * yk[1];

        u[0] = s[0] * matriz_jacobiana[0][0] + s[1] * matriz_jacobiana[1][0];
        u[1] = s[0] * matriz_jacobiana[0][1] + s[1] * matriz_jacobiana[1][1];

        producto_punto = u[0] * z[0] + u[1] * z[1];
        iteracion++;

        if (fabs(producto_punto) < tolerancia) {
            printf("Error Broyden: La matriz no es invertible. El método diverge.\n");
            return;
        }

        matriz_jacobiana[0][0] += (s[0] + z[0]) * u[0] / producto_punto;
        matriz_jacobiana[0][1] += (s[0] + z[0]) * u[1] / producto_punto;
        matriz_jacobiana[1][0] += (s[1] + z[1]) * u[0] / producto_punto;
        matriz_jacobiana[1][1] += (s[1] + z[1]) * u[1] / producto_punto;

        s[0] = -matriz_jacobiana[0][0] * F[0] - matriz_jacobiana[0][1] * F[1];
        s[1] = -matriz_jacobiana[1][0] * F[0] - matriz_jacobiana[1][1] * F[1];
        solucion[0] += s[0];
        solucion[1] += s[1];

        if (fabs(s[0]) < tolerancia && fabs(s[1]) < tolerancia) {
            printf("Método Broyden: Convergencia alcanzada en %d iteraciones\n", iteracion);
            return;
        }
    }
    printf("Error Broyden: El método no converge después de %d iteraciones.\n", iteracion);
}



void gradiente_conjugado_nolineal(double *solucion, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double)) {
    int iteracion = 0;
    double F[2], direccion[2], gradiente[2];
    double gradiente_anterior[2], direccion_anterior[2];
    double alpha, beta;

    while (iteraciones--) {
        F[0] = f1(solucion[0], solucion[1]);
        F[1] = f2(solucion[0], solucion[1]);
        gradiente[0] = F[0], gradiente[1] = F[1];

        if (iteracion == 0) {
            direccion[0] = -gradiente[0];
            direccion[1] = -gradiente[1];
        } else {
            beta = (gradiente[0] * gradiente[0] + gradiente[1] * gradiente[1]) /
                   (gradiente_anterior[0] * gradiente_anterior[0] + gradiente_anterior[1] * gradiente_anterior[1]);
            direccion[0] = -gradiente[0] + beta * direccion_anterior[0];
            direccion[1] = -gradiente[1] + beta * direccion_anterior[1];
        }
        alpha = 0.3;

        solucion[0] += alpha * direccion[0];
        solucion[1] += alpha * direccion[1];

        iteracion++;
        if (fabs(gradiente[0]) < tolerancia && fabs(gradiente[1]) < tolerancia) {
            printf("Gradiente Conjugado NoLineal: Convergencia alcanzada en %d iteraciones\n", iteracion);
            return;
        }
        gradiente_anterior[0] = gradiente[0], gradiente_anterior[1] = gradiente[1];
        direccion_anterior[0] = direccion[0], direccion_anterior[1] = direccion[1];
    }
    printf("Error Gradiente Conjugado NoLineal: No convergió después de %d iteraciones", iteracion);
}




double *euler(int N, double a, double b, double y0, double (*f)(double, double))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * sizeof(double));
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        y[i] = y[i - 1] + h * f(a, y[i - 1]);
        a += h;
    }
    return y;
}
double *heun(int N, double a, double b, double y0, double (*f)(double, double))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * sizeof(double));
    double yDummy, a_anterior = a;
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        a += h;
        yDummy = y[i - 1] + h * f(a_anterior, y[i - 1]);
        y[i] = y[i - 1] + (h / 2) * (f(a_anterior, y[i - 1]) + f(a, yDummy));
        a_anterior = a;
    }
    return y;
}
double derivada_parcial(double (*f)(double, double), double x, double y, double h, int variable)
{
    if (variable == 0)
    {
        return (f(x + h, y) - f(x - h, y)) / (2 * h);
    }
    if (variable == 1)
    {
        return (f(x, y + h) - f(x, y - h)) / (2 * h);
    }
    else
    {
        printf("Error en la variable\n");
        return 0;
    }
}
double *taylor_segundo_orden(int N, double a, double b, double y0, double (*f)(double, double))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * sizeof(double));
    double derivada_x, derivada_y;
    double h_2 = (h * h) / 2;
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        derivada_x = derivada_parcial(f, a, y[i - 1], h, 0);
        derivada_y = derivada_parcial(f, a, y[i - 1], h, 1);
        y[i] = y[i - 1] + h * f(a, y[i - 1]) + h_2 * (derivada_x + derivada_y * f(a, y[i - 1]));
        a += h;
    }

    return y;
}
double *RK4(int N, double a, double b, double y0, double (*f)(double, double))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * sizeof(double));
    double k1, k2, k3, k4;
    y[0] = y0;
    for (int i = 1; i <= N; i++)
    {
        k1 = h * f(a, y[i - 1]);
        k2 = h * f(a + h / 2, y[i - 1] + k1 / 2);
        k3 = h * f(a + h / 2, y[i - 1] + k2 / 2);
        k4 = h * f(a + h, y[i - 1] + k3);
        y[i] = y[i - 1] + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        a += h;
    }
    return y;
}


double *lotka_volterra(double t, double *xy)
{
    double *derivatives = malloc(2 * sizeof(double));
    derivatives[0] = 0.4 * xy[0] - 0.018 * xy[0] * xy[1];    
    derivatives[1] = -0.8 * xy[1] + 0.023 * xy[0] * xy[1];   
    return derivatives;
}

double *euler_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * dim * sizeof(double));
    double *dy;
    
      
    for (int j = 0; j < dim; j++)
    {
        y[j] = y0[j];
    }

    for (int i = 1; i <= N; i++)
    {
        dy = f(a, &y[(i - 1) * dim]);   
        for (int j = 0; j < dim; j++)
        {
            y[i * dim + j] = y[(i - 1) * dim + j] + h * dy[j];   
        }
        a += h;
        free(dy);   
    }
    return y;
}

double *heun_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * dim * sizeof(double));
    double *yDummy = malloc(dim * sizeof(double));
    double *f1, *f2;

      
    for (int j = 0; j < dim; j++)
    {
        y[j] = y0[j];
    }

    for (int i = 1; i <= N; i++)
    {
        f1 = f(a, &y[(i - 1) * dim]);
        for (int j = 0; j < dim; j++)
        {
            yDummy[j] = y[(i - 1) * dim + j] + h * f1[j];
        }
        f2 = f(a + h, yDummy);

        for (int j = 0; j < dim; j++)
        {
            y[i * dim + j] = y[(i - 1) * dim + j] + (h / 2) * (f1[j] + f2[j]);
        }
        a += h;
        free(f1), free(f2);
    }

    free(yDummy);
    return y;
}

double *RK4_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * dim * sizeof(double));
    double *k1, *k2, *k3, *k4, *yTmp;
    yTmp = malloc(dim * sizeof(double));

      
    for (int j = 0; j < dim; j++)
    {
        y[j] = y0[j];
    }

    for (int i = 1; i <= N; i++)
    {
        k1 = f(a, &y[(i - 1) * dim]);
        for (int j = 0; j < dim; j++)
        {
            yTmp[j] = y[(i - 1) * dim + j] + 0.5 * k1[j];
        }
        k2 = f(a + 0.5 * h, yTmp);

        for (int j = 0; j < dim; j++)
        {
            yTmp[j] = y[(i - 1) * dim + j] + 0.5 * k2[j];
        }
        k3 = f(a + 0.5 * h, yTmp);

        for (int j = 0; j < dim; j++)
        {
            yTmp[j] = y[(i - 1) * dim + j] + k3[j];
        }
        k4 = f(a + h, yTmp);

        for (int j = 0; j < dim; j++)
        {
            y[i * dim + j] = y[(i - 1) * dim + j] + (1.0 / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        }

        a += h;
        free(k1);
        free(k2);
        free(k3);
        free(k4);
    }
    free(yTmp);
    return y;
}

double *derivada_parcial_sistema(double *(*f)(double, double *), double x, double *y, double h, int variable)
{
    double yPlus[2], yMinus[2], *fPlus, *fMinus, *derivatives;
    derivatives = malloc(2 * sizeof(double));

    yPlus[0] = y[0];
    yPlus[1] = y[1];
    yMinus[0] = y[0];
    yMinus[1] = y[1];

    if (variable == 0)
    {
        yPlus[0] += h;
        yMinus[0] -= h;
    }
    else if (variable == 1)
    {
        yPlus[1] += h;
        yMinus[1] -= h;
    }
    else
    {
        printf("Error en la variable\n");
        return NULL;
    }

    fPlus = f(x, yPlus);
    fMinus = f(x, yMinus);

    derivatives[0] = (fPlus[0] - fMinus[0]) / (2 * h);   
    derivatives[1] = (fPlus[1] - fMinus[1]) / (2 * h);   

    free(fPlus);
    free(fMinus);

    return derivatives;
}

double *taylor_segundo_orden_sistema(int N, double a, double b, double *y0, int dim, double *(*f)(double, double *))
{
    double h = (b - a) / N;
    double *y = malloc((N + 1) * dim * sizeof(double));
    double *derivada_x, *derivada_y, *fy;
    double h_2 = (h * h) / 2;

      
    for (int j = 0; j < dim; j++)
    {
        y[j] = y0[j];
    }

    for (int i = 1; i <= N; i++)
    {
        fy = f(a, &y[(i - 1) * dim]);
        derivada_x = derivada_parcial_sistema(f, a, &y[(i - 1) * dim], h, 0);
        derivada_y = derivada_parcial_sistema(f, a, &y[(i - 1) * dim], h, 1);

        for (int j = 0; j < dim; j++)
        {
            y[i * dim + j] = y[(i - 1) * dim + j] + h * fy[j] + h_2 * (derivada_x[j] + derivada_y[j] * fy[j]);
        }

        a += h;
        free(fy);
        free(derivada_x);
        free(derivada_y);
    }
    return y;
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
