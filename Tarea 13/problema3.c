#include <stdio.h>
#include <stdlib.h>

double *lotka_volterra_derivatives(double t, double *xy) {
    double *derivatives = malloc(2 * sizeof(double));
    derivatives[0] = 0.4 * xy[0] - 0.018 * xy[0] * xy[1];  
    derivatives[1] = -0.8 * xy[1] + 0.023 * xy[0] * xy[1];  
    return derivatives;
}

double *euler(int N, double a, double b, double *y0, int dim, double *(*f)(double, double*)){
    double h = (b - a) / N;
    double *y = malloc((N+1) * dim * sizeof(double));
    double *dy;

     
    for(int j = 0; j < dim; j++){
        y[j] = y0[j];
    }

    for(int i = 1; i <= N; i++){
        dy = f(a, &y[(i-1) * dim]);  
        for(int j = 0; j < dim; j++){
            y[i*dim + j] = y[(i-1)*dim + j] + h * dy[j];  
        }
        a += h;
        free(dy);  
    }
    return y;
}

double *heun(int N, double a, double b, double *y0, int dim, double *(*f)(double, double*)){
    double h = (b - a) / N;
    double *y = malloc((N+1) * dim * sizeof(double));
    double *yDummy = malloc(dim * sizeof(double));
    double *f1, *f2;

     
    for(int j = 0; j < dim; j++){
        y[j] = y0[j];
    }

    for(int i = 1; i <= N; i++){
        f1 = f(a, &y[(i-1) * dim]);
        for(int j = 0; j < dim; j++){
            yDummy[j] = y[(i-1)*dim + j] + h * f1[j];
        }
        f2 = f(a + h, yDummy);

        for(int j = 0; j < dim; j++){
            y[i*dim + j] = y[(i-1)*dim + j] + (h/2)*(f1[j] + f2[j]);
        }
        a += h;
        free(f1), free(f2);
    }

    free(yDummy);

     
    printf("\nSolucion por Heun\n");
    for(int i = 0; i <= N; i++){
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, y[i*dim], y[i*dim + 1]);
    }

    return y;
}

double *RK4(int N, double a, double b, double *y0, int dim, double *(*f)(double, double*)) {
    double h = (b - a) / N;
    double *y = malloc((N+1) * dim * sizeof(double));
    double *k1, *k2, *k3, *k4, *yTmp;
    yTmp = malloc(dim * sizeof(double));

     
    for(int j = 0; j < dim; j++){
        y[j] = y0[j];
    }

    for(int i = 1; i <= N; i++){
        k1 = f(a, &y[(i-1) * dim]);
        for(int j = 0; j < dim; j++){
            yTmp[j] = y[(i-1)*dim + j] + 0.5 * k1[j];
        }
        k2 = f(a + 0.5 * h, yTmp);

        for(int j = 0; j < dim; j++){
            yTmp[j] = y[(i-1)*dim + j] + 0.5 * k2[j];
        }
        k3 = f(a + 0.5 * h, yTmp);

        for(int j = 0; j < dim; j++){
            yTmp[j] = y[(i-1)*dim + j] + k3[j];
        }
        k4 = f(a + h, yTmp);

        for(int j = 0; j < dim; j++){
            y[i*dim + j] = y[(i-1)*dim + j] + (1.0/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }

        a += h;
        free(k1); free(k2); free(k3); free(k4);
    }

    free(yTmp);

     
    printf("\nSolucion por Runge-Kutta de orden 4\n");
    for(int i = 0; i <= N; i++){
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, y[i*dim], y[i*dim + 1]);
    }

    return y;
}

double *derivada_parcial_x(double *(*f)(double, double*), double x, double *y, double h, int variable){
    double yPlus[2], yMinus[2], *fPlus, *fMinus, *derivatives;
    derivatives = malloc(2 * sizeof(double));
    
    yPlus[0] = y[0]; yPlus[1] = y[1];
    yMinus[0] = y[0]; yMinus[1] = y[1];

    if(variable == 0){
        yPlus[0] += h;
        yMinus[0] -= h;
    } else if (variable == 1){
        yPlus[1] += h;
        yMinus[1] -= h;
    } else {
        printf("Error en la variable\n");
        return NULL;
    }

    fPlus = f(x, yPlus);
    fMinus = f(x, yMinus);

    derivatives[0] = (fPlus[0] - fMinus[0]) / (2 * h);  
    derivatives[1] = (fPlus[1] - fMinus[1]) / (2 * h);  

    free(fPlus);
    free(fMinus);

    return derivatives;
}

 

double *taylor_segundo_orden(int N, double a, double b, double *y0, int dim, double *(*f)(double, double*)){
    double h = (b - a) / N;
    double *y = malloc((N+1) * dim * sizeof(double));
    double *derivada_x, *derivada_y, *fy;
    double h_2 = (h * h) / 2;
    
     
    for(int j = 0; j < dim; j++){
        y[j] = y0[j];
    }

    for(int i = 1; i <= N; i++){
        fy = f(a, &y[(i-1) * dim]);
        derivada_x = derivada_parcial_x(f, a, &y[(i-1) * dim], h, 0);
        derivada_y = derivada_parcial_x(f, a, &y[(i-1) * dim], h, 1);

        for(int j = 0; j < dim; j++){
            y[i*dim + j] = y[(i-1)*dim + j] + h*fy[j] + h_2*(derivada_x[j] + derivada_y[j]*fy[j]);
        }

        a += h;
        free(fy); free(derivada_x); free(derivada_y);
    }

     
    printf("\nSolucion por Taylor de segundo orden\n");
    for(int i = 0; i <= N; i++){
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, y[i*dim], y[i*dim + 1]);
    }

    return y;
}

int main() {
    int N = 1000;  
    double a = 0, b = 25;  
    double y0[] = {30, 4};  
    int dim = 2;  

    double *resultados = euler(N, a, b, y0, dim, lotka_volterra_derivatives);

     
    for(int i = 0; i <= N; i++) {
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }

     
    FILE *fp = fopen("lotka_volterra_euler.txt", "w");
    fprintf(fp, "t x y\n");
    for(int i = 0; i <= N; i++) {
        fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
    fclose(fp);

     
    resultados = heun(N, a, b, y0, dim, lotka_volterra_derivatives);
    
     
    for(int i = 0; i <= N; i++) {
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
     
    fp = fopen("lotka_volterra_heun.txt", "w");
    fprintf(fp, "t x y\n");
    for(int i = 0; i <= N; i++) {
        fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
    fclose(fp);

     
    N = 48;
    resultados = RK4(N, a, b, y0, dim, lotka_volterra_derivatives);

     
    for(int i = 0; i <= N; i++) {
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
     
    fp = fopen("lotka_volterra_RK4.txt", "w");
    fprintf(fp, "t x y\n");
    for(int i = 0; i <= N; i++) {
        fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
    fclose(fp);

     
    N = 1000;
    resultados = taylor_segundo_orden(N, a, b, y0, dim, lotka_volterra_derivatives);

     
    for(int i = 0; i <= N; i++) {
        printf("t = %f, x = %f, y = %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
     
    fp = fopen("lotka_volterra_taylor.txt", "w");
    fprintf(fp, "t x y\n");
    for(int i = 0; i <= N; i++) {
        fprintf(fp, "%f %f %f\n", a + i * (b - a) / N, resultados[i*dim], resultados[i*dim + 1]);
    }
    fclose(fp);



    free(resultados);  
    return 0;
}
