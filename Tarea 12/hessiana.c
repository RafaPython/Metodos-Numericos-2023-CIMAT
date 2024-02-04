#include <stdio.h>

void hessiana3(double **H, double x, double y, double z, double h,  double (*f)(double, double, double)) {
    H[0][0] = (f(x + h, y, z) - 2 * f(x, y, z) + f(x - h, y, z)) / (h * h);
    H[1][1] = (f(x, y + h, z) - 2 * f(x, y, z) + f(x, y - h, z)) / (h * h);
    H[2][2] = (f(x, y, z + h) - 2 * f(x, y, z) + f(x, y, z - h)) / (h * h);

    H[0][1] = H[1][0] = (f(x + h, y + h, z) - f(x + h, y - h, z) - f(x - h, y + h, z) + f(x - h, y - h, z)) / (4 * h * h);
    H[0][2] = H[2][0] = (f(x + h, y, z + h) - f(x + h, y, z - h) - f(x - h, y, z + h) + f(x - h, y, z - h)) / (4 * h * h);
    H[1][2] = H[2][1] = (f(x, y + h, z + h) - f(x, y + h, z - h) - f(x, y - h, z + h) + f(x, y - h, z - h)) / (4 * h * h);
}

void hessiana2(double **H, double x, double y, double h,  double (*f)(double, double)) {
    H[0][0] = (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / (h * h);
    H[1][1] = (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / (h * h);

    H[0][1] = H[1][0] = (f(x + h, y + h) - f(x + h, y - h) - f(x - h, y + h) + f(x - h, y - h)) / (4 * h * h);
}
