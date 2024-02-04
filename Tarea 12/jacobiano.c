#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void jacobiano3(double J[3][3], double x, double y, double z, double h, double (*f1)(double,double,double), double (*f2)(double,double,double), double (*f3)(double,double,double)){
    J[0][0] = (1/(12*h))*(-25*f1(x,y,z) + 48*f1(x+h,y,z) - 36*f1(x+2*h,y,z) + 16*f1(x+3*h,y,z) - 3*f1(x+4*h,y,z));
    J[0][1] = (1/(12*h))*(-25*f1(x,y,z) + 48*f1(x,y+h,z) - 36*f1(x,y+2*h,z) + 16*f1(x,y+3*h,z) - 3*f1(x,y+4*h,z));
    J[0][2] = (1/(12*h))*(-25*f1(x,y,z) + 48*f1(x,y,z+h) - 36*f1(x,y,z+2*h) + 16*f1(x,y,z+3*h) - 3*f1(x,y,z+4*h));
    J[1][1] = (1/(12*h))*(-25*f2(x,y,z) + 48*f2(x,y+h,z) - 36*f2(x,y+2*h,z) + 16*f2(x,y+3*h,z) - 3*f2(x,y+4*h,z));
    J[2][2] = (1/(12*h))*(-25*f3(x,y,z) + 48*f3(x,y,z+h) - 36*f3(x,y,z+2*h) + 16*f3(x,y,z+3*h) - 3*f3(x,y,z+4*h));
    J[1][2] = (1/(12*h))*(-25*f2(x,y,z) + 48*f2(x,y,z+h) - 36*f2(x,y,z+2*h) + 16*f2(x,y,z+3*h) - 3*f2(x,y,z+4*h));
    J[2][1] = (1/(12*h))*(-25*f3(x,y,z) + 48*f3(x,y+h,z) - 36*f3(x,y+2*h,z) + 16*f3(x,y+3*h,z) - 3*f3(x,y+4*h,z));
    J[1][0] = (1/(12*h))*(-25*f2(x,y,z) + 48*f2(x+h,y,z) - 36*f2(x+2*h,y,z) + 16*f2(x+3*h,y,z) - 3*f2(x+4*h,y,z));
    J[2][0] = (1/(12*h))*(-25*f3(x,y,z) + 48*f3(x+h,y,z) - 36*f3(x+2*h,y,z) + 16*f3(x+3*h,y,z) - 3*f3(x+4*h,y,z));   
}
void jacobiano2(double J[2][2], double x, double y, double h, double (*f1)(double,double), double (*f2)(double,double)){
    J[0][0] = (1/(12*h))*(-25*f1(x,y) + 48*f1(x+h,y) - 36*f1(x+2*h,y) + 16*f1(x+3*h,y) - 3*f1(x+4*h,y));
    J[0][1] = (1/(12*h))*(-25*f1(x,y) + 48*f1(x,y+h) - 36*f1(x,y+2*h) + 16*f1(x,y+3*h) - 3*f1(x,y+4*h));
    J[1][1] = (1/(12*h))*(-25*f2(x,y) + 48*f2(x,y+h) - 36*f2(x,y+2*h) + 16*f2(x,y+3*h) - 3*f2(x,y+4*h));
    J[1][0] = (1/(12*h))*(-25*f2(x,y) + 48*f2(x+h,y) - 36*f2(x+2*h,y) + 16*f2(x+3*h,y) - 3*f2(x+4*h,y));
}

void hessiana3(double **H, double x, double y, double z, double h, 
               double (*f)(double, double, double)) {
    // Segundas derivadas parciales
    H[0][0] = (f(x + h, y, z) - 2 * f(x, y, z) + f(x - h, y, z)) / (h * h);
    H[1][1] = (f(x, y + h, z) - 2 * f(x, y, z) + f(x, y - h, z)) / (h * h);
    H[2][2] = (f(x, y, z + h) - 2 * f(x, y, z) + f(x, y, z - h)) / (h * h);

    // Derivadas parciales mixtas
    H[0][1] = H[1][0] = (f(x + h, y + h, z) - f(x + h, y - h, z) - f(x - h, y + h, z) + f(x - h, y - h, z)) / (4 * h * h);
    H[0][2] = H[2][0] = (f(x + h, y, z + h) - f(x + h, y, z - h) - f(x - h, y, z + h) + f(x - h, y, z - h)) / (4 * h * h);
    H[1][2] = H[2][1] = (f(x, y + h, z + h) - f(x, y + h, z - h) - f(x, y - h, z + h) + f(x, y - h, z - h)) / (4 * h * h);
}

void hessiana2(double **H, double x, double y, double h, 
               double (*f)(double, double)) {
    // Segundas derivadas parciales
    H[0][0] = (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / (h * h);
    H[1][1] = (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / (h * h);

    // Derivada parcial mixta
    H[0][1] = H[1][0] = (f(x + h, y + h) - f(x + h, y - h) - f(x - h, y + h) + f(x - h, y - h)) / (4 * h * h);
}


double f1(double x, double y, double z) {
    return pow(x,4) + 3*pow(y,2)*x;
}

double f2(double x, double y, double z) {
    return 5*y*y - 2*x*y + 1;
}

double f3(double x, double y, double z) {
    return x + y + z;
}
double f12(double x, double y) {
    return pow(x,4) + 3*pow(y,2)*x;
}

double f22(double x, double y) {
    return 5*y*y - 2*x*y + 1;
}
double fh(double x, double y) {
    return x*x*y + y*y*x;
}
double fh1(double x, double y, double z) {
    return exp(-x)*sin(y*z);
}
int main() {
    // Parámetros de ejemplo
    double x = 1.0, y = 1.0, h = 0.001;

    // Crear matriz para el jacobiano
    // double J[3][3];

    // // Calcular el jacobiano
    // jacobiano3(J, x, y, z, h, f1, f2, f3);

    // // Imprimir el resultado
    // printf("Jacobiano 3 variables:\n");
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         printf("%.1f ", J[i][j]);
    //     }
    //     printf("\n");
    // }
    // double J2[2][2];

    // jacobiano2(J2, x, y, h, f12, f22);

    // printf("\nJacobiano 2 variables:\n");
    // for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 2; j++) {
    //         printf("%.1f ", J2[i][j]);
    //     }
    //     printf("\n");
    // }
    // calculamos para las hesianas
    double **H = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        H[i] = (double *) malloc(3 * sizeof(double));
    }
    double x1 = 0.0, y1 = 1.0, z1 = 3.14159265358979323846;
    hessiana3(H, x1, y1, z1, h, fh1);
    printf("\nHessiana 3 variables:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j <3; j++) {
            printf("%.4f ", H[i][j]);
        }
        printf("\n");
    }
    double **H2 = (double **) malloc(2 * sizeof(double *));
    for (int i = 0; i < 2; i++) {
        H2[i] = (double *) malloc(2 * sizeof(double));
    }
    hessiana2(H2, x, y, h, fh);
    printf("\nHessiana 2 variables:\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j <2; j++) {
            printf("%.1f ", H2[i][j]);
        }
        printf("\n");
    }
    // Liberar memoria
    for (int i = 0; i < 3; i++) {
        free(H[i]);
    }
    free(H);
    for (int i = 0; i < 2; i++) {
        free(H2[i]);
    }
    free(H2);

    // 

    return 0;
}