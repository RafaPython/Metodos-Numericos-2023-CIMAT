#include <stdio.h>
#include <math.h>

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
