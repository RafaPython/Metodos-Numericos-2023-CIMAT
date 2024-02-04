#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "solvers.h"


double **theta_metodo_calor(int N, int iteraciones, double theta, double a, double b, double t, double c, double (*f)(double)){
    
    double h = (b-a)/N;
    double lambda = (c*(t/N))/(h*h);

    double **A = generar_matriz(N);
    double **B = generar_matriz(N);
    double **C = generar_matriz(N);

    double **temp1 = generar_matriz(N);
    double **temp2 = generar_matriz(N);

    // generamos la matriz de calor 

    for (int i = 0; i < N; i++) {
        A[i][i] = 2.0;  // Diagonal principal
        B[i][i] = 1.0;  // Diagonal principal

        if (i < N - 1) {
            A[i][i + 1] = -1.0;  // Diagonal superior
            A[i + 1][i] = -1.0;  // Diagonal inferior
        }
    }

    for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
            if (i != j && i != j - 1 && i != j + 1) A[i][j] = 0.0;

    for (int i = 0; i < N; i++) C[0][i] = f(a + h*i);

    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) temp1[i][j] = theta*lambda*A[i][j];
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) temp1[i][j] += B[i][j];
    
    double dummy = theta - 1;
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) temp2[i][j] = dummy*lambda*A[i][j];
    for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) temp2[i][j] += B[i][j];   

    double tolerancia = 1e-5;
    for (int j = 1; j < N; j++) 
        C[j] = gauss_seidel_(N, temp1, multiplicar_matriz_vector_(N, temp2, C[j - 1]), tolerancia, iteraciones);

    // liberamos memoria 
    liberar_matriz(N,A);
    liberar_matriz(N, B);
    liberar_matriz(N, temp1);
    liberar_matriz(N, temp2);
    
    return C;
}

double intial_condition(double x){
    return 4*x - 4*x*x;
}

int main(){

    //  resolvemos el sistema 

    int N , iteraciones;
    double theta, a, b, t, c;
    
    printf("Ingrese el numero de puntos N: ");
    scanf("%d", &N);
    printf("Ingrese el numero de iteraciones: ");
    printf("Ingrese el valor de theta: ");
    scanf("%lf", &theta);
    printf("Ingrese el valor de a: ");
    scanf("%lf", &a);
    printf("Ingrese el valor de b: ");
    scanf("%lf", &b);
    printf("Ingrese el valor de t: ");
    scanf("%lf", &t);
    printf("Ingrese el valor de c: ");
    scanf("%lf", &c);


    double **solucion = theta_metodo_calor(N, iteraciones, theta, a, b, t, c, intial_condition);

    //  guardamos la solucion en un archivo de texto
    FILE *fp;
    fp = fopen("theta.txt", "w");

    printf("Como desea guardar la solucion? \n");
    printf("1. Guardar en un archivo de texto la matriz completa \n");
    printf("2. Guardar en un archivo de texto la matriz en formato : i,j, solucion[i][j] (esta opcion graficara tambien con gnuplot)\n");
    int opcion; 
    scanf("%d", &opcion);

    if (opcion == 1){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                fprintf(fp, "%.2f\t", solucion[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    else if (opcion == 2){
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) fprintf(fp, "%d %d %f\n", i, j, solucion[i][j]);
            fprintf(fp, "\n");  // Línea en blanco después de cada fila
        }
        fclose(fp);
    }
    else{
        printf("Opcion no valida; guardando como formato 1\n");
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                fprintf(fp, "%.2f\t", solucion[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }


    // liberar_matriz(N, solucion);

    printf("Se ha guardado la solucion en el archivo theta.txt para los siguientes datos: \n");
    printf("N = %d \n", N);
    printf("theta = %f \n", theta);
    printf("a = %f \n", a);
    printf("b = %f \n", b);
    printf("t = %f \n", t);
    printf("c = %f \n", c);
    printf("iteraciones = %d \n", iteraciones);
    printf("Tolerancia = %f \n", 1e-5);
    
    if (opcion == 2){
           // Crear el archivo de script de gnuplot
        FILE *gp;
        gp = fopen("plot.gp", "w");
        fprintf(gp, "set terminal qt font \"Arial\"\n");
        fprintf(gp, "set title \"Gráfico solución ecuación de calor\"\n");
        fprintf(gp, "set xlabel \"Eje X\"\n");
        fprintf(gp, "set ylabel \"Eje Y\"\n");
        fprintf(gp, "set zlabel \"Solución\"\n");
        fprintf(gp, "splot \"theta.txt\" using 1:2:3 with pm3d\n");
        fprintf(gp, "pause mouse close\n");  
        fclose(gp);

        // Ejecutar el script de gnuplot
        system("gnuplot plot.gp");
    }
    return 0;
}