/******************************************************************************
                               TAREA 8
                        M�TODOS NUM�RICOS
           Programa M�todo Rayleigh de Factorizaci�n.
           Autor: Enya Tovar Estrada.    Fecha: 06.OCT.23.
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Ya se implementaron los calloc
///////////////////////////////////////////////////////////////////////////////////
//Funci�n para crear las matrices:
double **CrearMatriz(int m, int n){
    double **matriz = (double **)calloc(m, sizeof(double *));
    if (matriz != NULL){
        for (int i = 0; i < m; i++){
            matriz[i] = (double *)calloc(n, sizeof(double));
            if (matriz[i] == NULL){
                printf("\nError al asignar memoria para matriz.\n");
                exit(1);
            }
        }
    }
    else {
        return NULL;
        printf("\nError al asignar memoria para matriz.\n");
        exit(1);
    }
    return matriz;
}
///////////////////////////////////////////////////////////////////////////////////////
//Funcion para crear vectores:
double *CrearVector(int n){
    double *vector = (double *)calloc(n, sizeof(double));
    if (vector == NULL){
        printf("\nError al asignar memoria para vector.\n");
    }
    return vector;
}
///////////////////////////////////////////////////////////////////////////////////////
//Liberar memoria de las matrices:
void LiberarMatriz(double **matriz, int m){
    for(int i=0; i<m; i++){
            free(matriz[i]);
        }
free(matriz);
}
///////////////////////////////////////////////////////////////////////////////////
//Leer matriz desde archivo
void leerMatrizDesdeArchivo(char *nombreArchivo, double ***matriz, int *m, int *n){
    FILE *archivo;
    int i,j;
    archivo = fopen(nombreArchivo,"r");
    if (archivo != NULL){
        fscanf(archivo, "%d %d",m,n);
        *matriz = CrearMatriz(*m,*n);
        for (i=0; i<*m; i++){
            for (j = 0; j < *n; j++){
                fscanf(archivo, "%lf", &(*matriz)[i][j]);
                }
            }
        fclose(archivo);
        }
else{
    printf("El archivo no se abrio\n");
    exit(1);//No se pudo abrir el archivo.
    }
}
//////////////////////////////////////////////////////////////////////////////////////
//Escribir el vector en archivo de texto
void EscribirVector(char *nombreArchivo, double *vector, int n){
    FILE *archivo;
    archivo = fopen(nombreArchivo, "w");
    int m=1;
    if (archivo != NULL) {
        fprintf(archivo, "%d %d\n",n,m);
        for (int i=0; i < n; i++){
            fprintf(archivo,"%.6lf\n",vector[i]);
        }
        fclose(archivo);
    } else {
        printf("No se pudo abrir el archivo para escritura.\n");
    }
}
/////////////////////////////////////////////////////////////////////////////////////////
//Funci�n para multiplicar Matriz Vector
void MultiplicarMatrizVector(double **A, int n, double *V0, double *V1){
    for (int i=0; i < n; i++) {
        double suma=0;
        for (int j=0; j < n; j++) {
            suma += A[i][j] * V0[j];
        }
        V1[i] = suma;
    }
}
///////////////////////////////////////////////////////////////////////////////////////
//Funci�n para normalizar vector
void NormalizarVector(double *V, int n){
    double suma=0;
    for (int i=0; i < n; i++) {
        suma += V[i]*V[i];
    }
    double magnitud = sqrt(suma);
    for (int i=0; i < n; i++) {
        V[i]/= magnitud;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
//Crear Matriz Identidad:
void MatrizIdentidad(int m, int n, double **matriz){
    if (n <= 0 || m <= 0){
        printf("No hay dimensiones para la I.\n");
        return;
    }
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                matriz[i][j] = 1;
            }
            else{
                matriz[i][j] = 0;
                }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////
//Producto punto vectores
double ProductoEscalar(double *x, double *y, int n){
    double resultado = 0;
    for (int i = 0; i < n; i++) {
        resultado += x[i]*y[i];
    }
    return resultado;
}
//////////////////////////////////////////////////////////////////////////////////////
// Funci�n para copiar un vector en otro:
void CopiarVector(double *v1, double *v2, int n){
    for(int i=0; i<n; i++){
        v2[i] = v1[i];
    }
}
//////////////////////////////////////////////////////////////////////////////////////////
//Factorizacion LU
void FactorizacionLU(double** A, double** L, double** U, int n){
    for (int i=0; i<n; i++){
        //U
        for(int k=i; k<n;k++){
            double suma=0;
            for (int j=0; j<i; j++){
                suma += L[i][j]*U[j][k];
            }
            U[i][k] = A[i][k] - suma;
        }
        //L
        for (int k=i; k<n; k++){
            if (i==k){
                L[i][i] = 1;
            }
            else{
                double suma = 0;
                for (int j=0; j<i; j++){
                    suma += L[k][j] * U[j][i];
                }
                L[k][i] = (A[k][i]-suma)/U[i][i];
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////
// Sustitucion de abajo hacia arriba y de arriba hacia abajo para resolver LUx = b
void ResolverSistemaLinealLU(double **L, double **U, double *b, double *x, int n){
    double *y = (double*)malloc(n*sizeof(double));

    //Ly=b
    for (int i=0; i<n; i++){
        double suma = 0;
        for (int j=0; j<i; j++){
            suma += L[i][j]*y[j];
        }
        y[i]=(b[i]-suma)/L[i][i];
    }
    //Ux=y
    for (int i=n-1; i>= 0; i--){
        double suma=0;
        for (int j = i+1; j<n; j++){
            suma += U[i][j]*x[j];
        }
        x[i]=(y[i]-suma)/U[i][i];
    }
    free(y);
}
/////////////////////////////////////////////////////////////////////////////////////////
//Aqu� empieza el m�todo de Rayleigh:
double *Rayleigh(double **A, double *v0, int m, int n, double sigma){
    double TOL=1E-8;
    int MaxIter=10000;
    int k=0;
    v0[0]=1;
    double **B = CrearMatriz(m,n);
    double **I = CrearMatriz(m,n);
    double *v1 = CrearVector(m);
    double **L = CrearMatriz(m,n);
    double **U = CrearMatriz(m,n);
    double *Mult =CrearVector(m);
    double *Diferencia = CrearVector(m);
    //Hacemos la matriz identidad
    MatrizIdentidad(m,n,I);
    for(k=0;k<MaxIter;k++){
            for(int i=0;i<m;i++){
                for(int j=0; j<n;j++){
                    B[i][j]= A[i][j] - sigma*I[i][j];
                }
            }
            //Resolver Bv1=v0
            FactorizacionLU(B,L,U,m);

            // revisamos que B = LU y que L y U sean triangulares
            double **aux = CrearMatriz(m,n);
            MultiplicarMatriz(L,U,m,aux);
            printf("Matriz B:\n");
            ImprimirMatriz(B,m);
            printf("\n");
            exit(0);

            ResolverSistemaLinealLU(L, U, v0, v1, m);
            NormalizarVector(v1,m);

            for(int i=0;i<m;i++){
                    Diferencia[i] = v1[i]- v0[i];
                    }
            double norma = sqrt(ProductoEscalar(Diferencia,Diferencia,m));
            if (norma < TOL){
                printf("\nConverge a %d iteraciones.\n",k);
                LiberarMatriz(L,m);
                LiberarMatriz(U,m);
                LiberarMatriz(B,m);
                LiberarMatriz(I,m);
                free(Diferencia);
                free(Mult);
                return(v1);
                break;
                }
            CopiarVector(v1,v0,m);
            MultiplicarMatrizVector(B,m,v0,Mult);
            sigma = ProductoEscalar(v0,Mult,m);
            }
printf("No converge a %d iteraciones.",k);
return(v1);
}
///////////////////////////////////////////////////////////////////////////////////////////
//Imprimir Matriz:
void ImprimirMatriz(double **A,int n){
    printf("Tu Matriz A es:\n");
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            printf("%.4lf ", A[i][j]);
        }
        printf("\n");
    }
}
///////////////////////////////////////////////////////////////////////////////////////////
//Imprimir vector:
void ImprimirVector(double *v, int n){
    for(int i=0; i<n; i++){
        printf("%.4lf \n",v[i]);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//PROGRAMA PRINCIPAL DONDE SOLO LLAMAMOS LAS FUNCIONES (ya mas chiquito)
int main(){
    printf("METODO DE RAYLEIGH PARA OBTENER LOS EIGENVALORES.\n");
    int m, n;
    double **A;
    double *v0,*v1;
    double sigma=1;
    v0 = CrearVector(m);
    v1 = CrearVector(m);
    leerMatrizDesdeArchivo("MatrixEjemplo.txt",&A,&m,&n); //Archivo entrada
    ImprimirMatriz(A,n);
    //Llamamos el m�todo:
    v1 = Rayleigh(A,v0,m,n,sigma);
    //Imprimimos soluci�n
    printf("Vector solucion:\n");
    ImprimirVector(v1,m);
    LiberarMatriz(A, m);
    free(v1);
    free(v0);
    return (0);
}
