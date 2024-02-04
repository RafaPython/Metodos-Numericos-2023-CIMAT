#include <stdio.h>
#include <stdlib.h>

void newton_nolineal(int N, double **J, double *F, double *x0, double tolerancia, double h, int iteraciones, double (*f1)(double,double), double (*f2)(double,double)){
    int i_convergencia = 0;
    double *x_anterior = (double *) malloc(N * sizeof(double));
    double dummy;
    while(iteraciones--){
        F[0] = f1(x0[0], x0[1]);
        F[1] = f2(x0[0], x0[1]);
        jacobiano2(J, x0[0], x0[1], h, f1, f2);

        dummy = -J[1][0] / J[0][0];
        J[1][1] -= dummy * J[0][1];
        F[1] -= dummy * F[0];

        x_anterior[1] = F[1] / J[1][1];
        x_anterior[0] = (F[0] - J[0][1] * x_anterior[1]) / J[0][0];

        x0[0] += x_anterior[0];
        x0[1] += x_anterior[1];

        if(fabs(x_anterior[0]) < tolerancia && fabs(x_anterior[1]) < tolerancia){
            printf("Newton no lineal: convergencia alcanzada en %d iteraciones\n", i_convergencia);
            free(x_anterior);
            return;
        }
        i_convergencia++;
    }
    printf("Newton no lineal: no converge en %d iteraciones\n", i_convergencia);
    free(x_anterior);
}