#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #include "solvers.h"


#include <stdlib.h>

// double spline_condicionado(int N, double *x, double *y, double *a, double *b, double *c, double *d, double punto) {
//     // Calcular los coeficientes 'h' y 'alpha'
//     double *h = (double *)malloc((N - 1) * sizeof(double));
//     double *alpha = (double *)malloc((N - 1) * sizeof(double));

//     for (int i = 0; i < N - 1; i++) {
//         h[i] = x[i + 1] - x[i];
//         alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / (i == 0 ? h[i] : h[i - 1])) * (y[i] - y[i - 1]);
//     }

//     // Calcular los coeficientes l, u y w
//     double *l = (double *)malloc(N * sizeof(double));
//     double *u = (double *)malloc(N * sizeof(double));
//     double *w = (double *)malloc(N * sizeof(double));

//     l[0] = 1.0; 
//     u[0] = 0.0;
//     w[0] = 0.0; 

//     for (int i = 1; i < N - 1; i++) {
//         l[i] = 2 * (h[i - 1] + h[i]) - h[i - 1] * u[i - 1];
//         u[i] = h[i] / l[i];
//         w[i] = (alpha[i] - alpha[i - 1]) - h[i - 1] * w[i - 1];
//         w[i] /= l[i];
//     }

//     l[N - 1] = 1.0; 
//     u[N - 1] = 0.0;
//     w[N - 1] = 0.0; 

//     // Calcular los coeficientes c, b, d y a
//     for (int i = N - 2; i >= 0; i--) {
//         c[i] = w[i] - u[i] * c[i + 1];
//         b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3.0;
//         d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
//         a[i] = y[i];
//     }

//     // Encontrar el intervalo correspondiente al punto
//     int intervalo = 0;
//     for (int i = 0; i < N - 1; i++) {
//         if (x[i] <= punto && punto <= x[i + 1]) {
//             intervalo = i;
//             break;
//         }
//     }

//     // Calcular el valor del spline en el punto
//     double dx = punto - x[intervalo];
//     double resultado = a[intervalo] + b[intervalo] * dx + c[intervalo] * dx * dx + d[intervalo] * dx * dx * dx;

//     // Liberar la memoria
//     free(h);
//     free(alpha);
//     free(l);
//     free(u);
//     free(w);

//     return resultado;
// }

#include <stdio.h>
#include <stdlib.h>

double spline_condicionado(int N, double *x, double *y, double *a, double *b, double *c, double *d, double *puntos_condicionados, double *valores_condicionados, double punto) {
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
    free(h);
    free(alpha);
    free(l);
    free(u);
    free(w);

    return resultado;
}


void imprimir_tabla(int N, double *x, double *a, double *b, double *c, double *d) {
    printf("i\tx[i]\t\ta[i]\t\tb[i]\t\tc[i]\t\td[i]\n");
    for (int i = 0; i < N; i++) {
        printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, x[i], a[i], b[i], c[i], d[i]);
    }
}

int main(){
    // exercise 28 of interpolation methods in the book Burden-Faires
    // exercise dog of the book Burden-Faires


    // we create the data for curve 1
    int N1 = 9;
    double x1[] = {1,2,5,6,7,8,10, 13, 17};
    double y1[] = {3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5};
    // puntos condicionados
    double puntos_condicionados1[] = {1, 17};
    // valores condicionados
    double valores_condicionados1[] = {1.0,-0.67};

    // data for curve 2
    int N2 = 7;
    double x2[] = {17, 20, 23, 24, 25, 27, 27.7};
    double y2[] = {4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1};
    // puntos condicionados
    double puntos_condicionados2[] = {17, 27.7};
    // valores condicionados
    double valores_condicionados2[] = {3.0, 0,0,0,0,0,4.0};

    // we create the data for curve 3
    int N3 = 4;
    double x3[] = {27.7, 28, 29, 30};
    double y3[] = {4.1, 4.3, 4.1, 3.0};
    // puntos condicionados
    double puntos_condicionados3[] = {27.7, 30};
    // valores condicionados
    double valores_condicionados3[] = {0.33, 0,0,-1.5};

    int punto1 = 3;
    int punto2 = 18;
    int punto3 = 28.5;

    double *a1 = (double *)malloc(N1 * sizeof(double));
    double *b1 = (double *)malloc(N1 * sizeof(double));
    double *c1 = (double *)malloc(N1 * sizeof(double));
    double *d1 = (double *)malloc(N1 * sizeof(double));

    double *a2 = (double *)malloc(N2 * sizeof(double));
    double *b2 = (double *)malloc(N2 * sizeof(double));
    double *c2 = (double *)malloc(N2 * sizeof(double));
    double *d2 = (double *)malloc(N2 * sizeof(double));

    double *a3 = (double *)malloc(N3 * sizeof(double));
    double *b3 = (double *)malloc(N3 * sizeof(double));
    double *c3 = (double *)malloc(N3 * sizeof(double));
    double *d3 = (double *)malloc(N3 * sizeof(double));


    // interpolamos para la curva 1
    double P1 = spline_condicionado(N1, x1, y1, a1, b1, c1, d1, puntos_condicionados1, valores_condicionados1, punto1);
    printf("Curva 1:\n");
    // imprimimos los coeficientes
    imprimir_tabla(N1, x1, a1, b1, c1, d1);

    // interpolamos para la curva 2
    printf("\n\n\n");
    double P2 = spline_condicionado(N2, x2, y2, a2, b2, c2, d2, puntos_condicionados2, valores_condicionados2, punto2);
    printf("Curva 2:\n");
    // imprimimos los coeficientes
    imprimir_tabla(N2, x2, a2, b2, c2, d2);

    // interpolamos para la curva 3
    printf("\n\n\n");
    double P3 =  spline_condicionado(N3, x3, y3, a3, b3, c3, d3, puntos_condicionados3, valores_condicionados3, punto3);
    printf("Curva 3:\n");
    // imprimimos los coeficientes
    imprimir_tabla(N3, x3, a3, b3, c3, d3);
    

    free(a1);
    free(b1);
    free(c1);
    free(d1);

    free(a2);
    free(b2);
    free(c2);
    free(d2);

    free(a3);
    free(b3);
    free(c3);
    free(d3);






    return 0;
}