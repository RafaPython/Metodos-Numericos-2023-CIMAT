#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void f(double x, double *y, double *salida) {
    salida[0] = y[1];
    salida[1] = 2 * y[1] - y[0] + sin(x);
}
double *RK4_sistema(int N, double a, double b, double *y0, int dim, void (*f)(double, double *, double *)) {
    double h = (b - a) / N;
    double *y = (double *)malloc((N + 1) * dim * sizeof(double));
    double k1[dim], k2[dim], k3[dim], k4[dim], yTmp[dim];
    for (int j = 0; j < dim; j++) y[j] = y0[j];
    for (int i = 1; i <= N; i++) {
        f(a, &y[(i - 1) * dim], k1);
        for (int j = 0; j < dim; j++) yTmp[j] = y[(i - 1) * dim + j] + 0.5 * h * k1[j];
        f(a + 0.5 * h, yTmp, k2);
        for (int j = 0; j < dim; j++) yTmp[j] = y[(i - 1) * dim + j] + 0.5 * h * k2[j];
        f(a + 0.5 * h, yTmp, k3);
        for (int j = 0; j < dim; j++) yTmp[j] = y[(i - 1) * dim + j] + h * k3[j];
        f(a + h, yTmp, k4);
        for (int j = 0; j < dim; j++) y[i * dim + j] = y[(i - 1) * dim + j] + (h / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        a += h;
    }
    return y;
}

int main(){
    double a = 0.0, b = 2.0;
    int N = 500;
    double *y0 = (double *)malloc(2 * sizeof(double));
    y0[0] = -1.0;
    y0[1] = 1.0;
    double *y = RK4_sistema(N, a, b, y0, 2, f);

    //  guardamos el resultado de x y en un archivo
    FILE *fp = fopen("problema4.txt", "w");
    for (int i = 0; i <= N; i++) 
        fprintf(fp, "%lf %lf\n", a + i * (b - a) / N, y[i * 2]);
    fclose(fp);
    free(y);
    free(y0);

    return 0;
}
