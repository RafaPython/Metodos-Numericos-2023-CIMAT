#include <stdio.h>

void Ux_b(int N, double U[N][N], double b[N] ,double x[N]);

int main(){
    const int N = 5;
    
    double U[5][5] = {
        {9.00000, 7.00000, 2.00000, 3.00000, 2.00000},
        {0.00000, 9.22222, 8.77778, 9.66667, 8.77778},
        {0.00000, 0.00000, 2.93976, -5.93976, 1.93976},
        {0.00000, 0.00000, 0.00000, -5.22951, 1.59016},
        {0.00000, 0.00000, 0.00000, 0.00000, 5.69749}
    };

    double b[5] = {51.0000, 127.3333, -5.2410, -12.9672, 28.4875};


    printf("d)\n\n");
    printf("Matriz Original :\n");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%f ", U[i][j]);
        printf("\n");
    }
    printf("Vector b:\n");
    for (int i = 0; i < N; i++) printf("%f\n", b[i]);

    double x[N];
    Ux_b(N, U, b, x);
    printf("Soluciones:\n");
    for (int i = 0; i < N; i++) printf("x[%d] = %f\n", i, x[i]);
    return 0;
}
void Ux_b(int N, double U[N][N], double b[N] ,double x[N]){
    double acumulador;
    for(int i = N - 1; i >= 0; i--){
        acumulador = 0;
        for (int k =  N; k > i; k--) acumulador += U[i][k]*x[k];
        x[i] = (b[i] - acumulador) / U[i][i];
    }  
}

