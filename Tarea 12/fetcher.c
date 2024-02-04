#include <stdio.h>
#include <math.h>

void gradiente_conjugado_nolineal(double *solucion, double tolerancia, int iteraciones, double (*f1)(double, double), double (*f2)(double, double)) {
    int iteracion = 0;
    double F[2], direccion[2], gradiente[2];
    double gradiente_anterior[2], direccion_anterior[2];
    double alpha, beta;

    while (iteraciones--) {
        F[0] = f1(solucion[0], solucion[1]);
        F[1] = f2(solucion[0], solucion[1]);
        gradiente[0] = F[0], gradiente[1] = F[1];

        if (iteracion == 0) {
            direccion[0] = -gradiente[0];
            direccion[1] = -gradiente[1];
        } else {
            beta = (gradiente[0] * gradiente[0] + gradiente[1] * gradiente[1]) /
                   (gradiente_anterior[0] * gradiente_anterior[0] + gradiente_anterior[1] * gradiente_anterior[1]);
            direccion[0] = -gradiente[0] + beta * direccion_anterior[0];
            direccion[1] = -gradiente[1] + beta * direccion_anterior[1];
        }
        alpha = 0.3;

        solucion[0] += alpha * direccion[0];
        solucion[1] += alpha * direccion[1];

        iteracion++;
        if (fabs(gradiente[0]) < tolerancia && fabs(gradiente[1]) < tolerancia) {
            printf("Gradiente Conjugado NoLineal: Convergencia alcanzada en %d iteraciones\n", iteracion);
            return;
        }
        gradiente_anterior[0] = gradiente[0], gradiente_anterior[1] = gradiente[1];
        direccion_anterior[0] = direccion[0], direccion_anterior[1] = direccion[1];
    }
    printf("Error Gradiente Conjugado NoLineal: No convergió después de %d iteraciones", iteracion);
}
