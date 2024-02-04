#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #include "solvers.h"

#include <stdio.h>
#include <stdlib.h>

double spline_natural(int N, double *x, double *y, double *a, double *b, double *c, double *d, double punto) {
    // double *a = (double *)malloc(N * sizeof(double));
    // double *b = (double *)malloc(N * sizeof(double));
    // double *c = (double *)malloc(N * sizeof(double));
    // double *d = (double *)malloc(N * sizeof(double));

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

    // free(a);
    // free(b);
    // free(c);
    // free(d);
    free(h);
    free(alpha);
    free(l);
    free(u);
    free(w);

    return spline;
}



int main(){
    // usamos el ejemplo del pato del libro Burden
    int N = 21;
    double x[] = {0.9, 1.3, 1.9, 2.1, 2.6, 3.0, 3.9, 4.4, 4.7, 5.0, 6.0, 7.0, 8.0, 9.2, 10.5, 11.3, 11.6, 12.0, 12.6, 13.0, 13.3};
    double y[] = {1.3, 1.5, 1.85, 2.1, 2.6, 2.7, 2.4, 2.15, 2.05, 2.1, 2.25, 2.3, 2.25, 1.95, 1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25};

    double *a = (double *)malloc(N * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *c = (double *)malloc(N * sizeof(double));
    double *d = (double *)malloc(N * sizeof(double));

    double puntos = 5.5;
    double P = spline_natural(N, x, y, a, b, c, d, puntos);
    printf("El valor de la interpolación de Spline Natural es: %lf\n", P);

    // imprimimos los coeficientes
    printf("j\txj\taj\tbj\tcj\tdj\n");
    for(int i = 0; i < N; i++) printf("%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n", i, x[i], a[i], b[i], c[i], d[i]);
    free(a);
    free(b);
    free(c);
    free(d);
    return 0;
}