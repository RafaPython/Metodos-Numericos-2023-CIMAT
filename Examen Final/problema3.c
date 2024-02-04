#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double interpolador_lagrange(int N, double *x, double *y, double punto){
    double L, P = 0.0;
    for (int i = 0; i < N; i++){
        L = 1;
        for (int j = 0; j < N; j++) if (i != j) L *= (punto - x[j])/(x[i] - x[j]);
        P += y[i] * L;
    }
    return P;
}

void interpolacion_cuadrados_exponencial(int N, double *x, double *y, double *alfa, double *beta, double regulador){
    double x1, suma_x = 0, suma_y = 0, suma_xy = 0, suma_x2 = 0;

    for (int i = 0; i < N; i++) {
        x1 = exp(-3 * x[i]);
        suma_x += x1;
        suma_y += y[i];
        suma_xy += x1 * y[i];
        suma_x2 += x1 * x1;
    }

    suma_x2 += regulador * N;
    double xx = suma_x2 - (suma_x * suma_x) / (double)N;
    double xy = suma_xy - (suma_x * suma_y) / (double)N;
    *beta = xy / xx;
    *alfa = (suma_y - (*beta) * suma_x) / (double) N;
}

int main(){
    int N = 4;
    double t[4] = {0.0, 3.0, 5.0, 7.0};
    double c[4] = {1.0, 20.0, 22.0, 23.0};
    //  calculamos con interpolación de lagrange t = 4
    double t0 = 4.0;
    double c0 = interpolador_lagrange(4, t, c, t0);
    printf("El valor de c(4) es: %lf\n", c0);
    // imprimimos el polinomio de interpolación
    printf("El polinomio de interpolación es:\n");
    printf("P(x) = %lf + %lf(x - %lf) + %lf(x - %lf)(x - %lf) + %lf(x - %lf)(x - %lf)(x - %lf)\n", c[0], (c[1] - c[0])/(t[1] - t[0]), t[0], (c[2] - c[1])/(t[2] - t[1]), t[1], t[0], (c[3] - c[2])/(t[3] - t[2]), t[2], t[1], t[0]);
    printf("\n\n");
    double w[] = {0, 0}; 
    double regulador = 0.0001; 
    for (int i = 0; i < N; i++) c[i] = 1 / c[i];
    interpolacion_cuadrados_exponencial(N, t, c, &w[0], &w[1], regulador);
    printf("Parametros encontrados: alpha = %f, beta = %f\n", w[0], w[1]);

    return 0;
}