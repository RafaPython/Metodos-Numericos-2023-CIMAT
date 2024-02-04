#include <stdio.h>

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