#include <stdio.h>

void Lx_b(int N, double L[N][N], double b[N] ,double x[N]);

int main(){
    const int N = 5;

    double L[5][5] = {
        {1.00000, 0.00000, 0.00000, 0.00000, 0.00000},
        {0.11111, 1.00000, 0.00000, 0.00000, 0.00000},
        {0.33333, 0.61446, 1.00000, 0.00000, 0.00000},
        {0.33333, 0.28916, -0.40984, 1.00000, 0.00000},
        {0.88889, 0.19277, -0.84016, 0.29075, 1.00000}
    };

    double b[5] = {51, 133, 90, 43, 99};





    printf("c)\n\n");
    printf("Matriz Original :\n");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%f ", L[i][j]);
        printf("\n");
    }
    printf("Vector b:\n");
    for (int i = 0; i < N; i++) printf("%f\n", b[i]);

    double x[N];
    Lx_b(N, L, b, x);
    printf("Soluciones:\n");
    for (int i = 0; i < N; i++) printf("x[%d] = %f\n", i, x[i]);
    return 0;
}
void Lx_b(int N, double L[N][N], double b[N] ,double x[N]){
    double acumulador;
    for(int i = 0; i < N; i++){
        acumulador = 0;
        for (int k = 0; k <= i; k++) acumulador += L[i][k]*x[k];
        x[i] = (b[i] - acumulador) / L[i][i];
    }  
}

