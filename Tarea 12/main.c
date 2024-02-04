#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "solvers.h"

double f1(double x, double y, double z)
{
    return (x + y - 3);
}
double f2(double x, double y, double z)
{
    return (x * x + y * y - 9);
}
double f3(double x, double y, double z)
{
    return (4 * y * y - 2 * z);
}

int main()
{
    int opcion;
    double x, y, z, h;
    double resultado;

    printf("Seleccione una opción:\n");
    printf("1. Calcular derivada hacia adelante\n");
    printf("2. Calcular derivada hacia atrás\n");
    printf("3. Calcular derivada centrada\n");
    printf("4. Calcular derivada de tres puntos en el extremo\n");
    printf("5. Calcular derivada de tres puntos en el punto medio\n");
    printf("6. Calcular derivada de cinco puntos en el punto medio\n");
    printf("7. Calcular derivada de cinco puntos en el extremo\n");
    printf("8. Calcular segunda derivada\n");
    printf("9. Calcular derivada hacia adelante con datos\n");
    printf("10. Calcular derivada hacia atrás con datos\n");
    printf("11. Calcular derivada centrada con datos\n");
    printf("12. Calcular derivada de tres puntos en el extremo con datos\n");
    printf("13. Calcular derivada de tres puntos en el punto medio con datos\n");
    printf("14. Calcular derivada de cinco puntos en el punto medio con datos\n");
    printf("15. Calcular derivada de cinco puntos en el extremo con datos\n");
    printf("16. Calcular segunda derivada con datos\n");
    printf("17. Calcular Jacobiano 3x3\n");
    printf("18. Calcular Jacobiano 2x2\n");
    printf("19. Calcular Hessiana 3x3\n");
    printf("20. Calcular Hessiana 2x2\n");
    printf("21. Calcular Punto Fijo No Lineal\n");
    printf("22. Calcular Método de Newton No Lineal\n");
    printf("23. Calcular Método de Broyden\n");
    printf("24. Calcular Gradiente Conjugado No Lineal\n");
    printf("25. Salir\n");

    scanf("%d", &opcion);

    if (opcion < 1 || opcion > 25)
    {
        printf("Opción no válida.\n");
        return 0;
    }

    switch (opcion)
    {
    case 1:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_hacia_adelante(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 2:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_hacia_atras(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 3:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_centrada(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 4:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_tres_puntos_endpoint(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 5:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_tres_puntos_midpoint(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 6:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_cinco_puntos_midpoint(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 7:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = derivada_cinco_puntos_endpoint(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 8:
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h);
        resultado = segunda_derivada(f1, x, h);
        printf("Resultado: %lf\n", resultado);
        break;
    case 9:
        // Calcular derivada hacia adelante con datos
        double f1_data, f2_data, h_data;
        printf("Ingrese el valor de f(x): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_hacia_adelante_datos(f1_data, f2_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 10:
        // Calcular derivada hacia atrás con datos
        printf("Ingrese el valor de f(x-h): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_hacia_atras_datos(f1_data, f2_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 11:
        // Calcular derivada centrada con datos
        printf("Ingrese el valor de f(x-h): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_centrada_datos(f1_data, f2_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 12:
        // Calcular derivada de tres puntos en el extremo con datos
        double f3_data;
        printf("Ingrese el valor de f(x): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de f(x+2h): ");
        scanf("%lf", &f3_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_tres_puntos_endpoint_datos(f1_data, f2_data, f3_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 13:
        // Calcular derivada de tres puntos en el punto medio con datos
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x-h): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_tres_puntos_midpoint_datos(f1_data, f2_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 14:
        // Calcular derivada de cinco puntos en el punto medio con datos
        double f4_data;
        printf("Ingrese el valor de f(x-2h): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x-h): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f3_data);
        printf("Ingrese el valor de f(x+2h): ");
        scanf("%lf", &f4_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_cinco_puntos_midpoint_datos(f1_data, f2_data, f3_data, f4_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 15:
        // Calcular derivada de cinco puntos en el extremo con datos
        double f5_data;
        printf("Ingrese el valor de f(x): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de f(x+2h): ");
        scanf("%lf", &f3_data);
        printf("Ingrese el valor de f(x+3h): ");
        scanf("%lf", &f4_data);
        printf("Ingrese el valor de f(x+4h): ");
        scanf("%lf", &f5_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = derivada_cinco_puntos_endpoint_datos(f1_data, f2_data, f3_data, f4_data, f5_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 16:
        // Calcular segunda derivada con datos
        printf("Ingrese el valor de f(x-h): ");
        scanf("%lf", &f1_data);
        printf("Ingrese el valor de f(x): ");
        scanf("%lf", &f2_data);
        printf("Ingrese el valor de f(x+h): ");
        scanf("%lf", &f3_data);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_data);
        resultado = segunda_derivada_midpoint_datos(f1_data, f2_data, f3_data, h_data);
        printf("Resultado: %lf\n", resultado);
        break;
    case 17:
        // Calcular Jacobiano 3x3
        double J3x3[3][3];
        double x_jacobiano, y_jacobiano, z_jacobiano, h_jacobiano;
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x_jacobiano);
        printf("Ingrese el valor de y: ");
        scanf("%lf", &y_jacobiano);
        printf("Ingrese el valor de z: ");
        scanf("%lf", &z_jacobiano);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_jacobiano);
        jacobiano3((double **)J3x3, x_jacobiano, y_jacobiano, z_jacobiano, h_jacobiano, f1, f2, f3);
        printf("Jacobiano 3x3:\n");
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                printf("%.6lf ", J3x3[i][j]);
            }
            printf("\n");
        }
        break;
    case 18:
        // Calcular Jacobiano 2x2
        double J2x2[2][2];
        double x_jacobiano_2x2, y_jacobiano_2x2, h_jacobiano_2x2;
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x_jacobiano_2x2);
        printf("Ingrese el valor de y: ");
        scanf("%lf", &y_jacobiano_2x2);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_jacobiano_2x2);
        jacobiano2((double **)J2x2, x_jacobiano_2x2, y_jacobiano_2x2, h_jacobiano_2x2, f1, f2);
        printf("Jacobiano 2x2:\n");
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                printf("%.6lf ", J2x2[i][j]);
            }
            printf("\n");
        }
        break;
    case 19:
        // Calcular Hessiana 3x3
        double H3x3[3][3];
        double x_hessiana, y_hessiana, z_hessiana, h_hessiana;
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x_hessiana);
        printf("Ingrese el valor de y: ");
        scanf("%lf", &y_hessiana);
        printf("Ingrese el valor de z: ");
        scanf("%lf", &z_hessiana);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_hessiana);
        hessiana3((double **)H3x3, x_hessiana, y_hessiana, z_hessiana, h_hessiana, f1, f2, f3);
        printf("Hessiana 3x3:\n");
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                printf("%.6lf ", H3x3[i][j]);
            }
            printf("\n");
        }
        break;
    case 20:
        // Calcular Hessiana 2x2
        double H2x2[2][2];
        double x_hessiana_2x2, y_hessiana_2x2, h_hessiana_2x2;
        printf("Ingrese el valor de x: ");
        scanf("%lf", &x_hessiana_2x2);
        printf("Ingrese el valor de y: ");
        scanf("%lf", &y_hessiana_2x2);
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_hessiana_2x2);
        hessiana2((double **)H2x2, x_hessiana_2x2, y_hessiana_2x2, h_hessiana_2x2, f1, f2);
        printf("Hessiana 2x2:\n");
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                printf("%.6lf ", H2x2[i][j]);
            }
            printf("\n");
        }
        break;
    case 21:
        // Calcular Punto Fijo No Lineal
        int N_punto_fijo;
        printf("Ingrese la cantidad de iteraciones: ");
        scanf("%d", &N_punto_fijo);
        double x0_punto_fijo[2];
        printf("Ingrese el valor inicial de x0: ");
        scanf("%lf", &x0_punto_fijo[0]);
        printf("Ingrese el valor inicial de x1: ");
        scanf("%lf", &x0_punto_fijo[1]);
        double tolerancia_punto_fijo;
        printf("Ingrese la tolerancia: ");
        scanf("%lf", &tolerancia_punto_fijo);
        punto_fijo_nolineal(2, x0_punto_fijo, tolerancia_punto_fijo, N_punto_fijo, f1, f2);
        printf("El valor de x1 es: %lf\n", x0_punto_fijo[0]);
        printf("El valor de x2 es: %lf\n", x0_punto_fijo[1]);
        break;
    case 22:
        // Calcular Método de Newton No Lineal
        int N_newton;
        printf("Ingrese la cantidad de iteraciones: ");
        scanf("%d", &N_newton);
        double x0_newton[2];
        printf("Ingrese el valor inicial de x0: ");
        scanf("%lf", &x0_newton[0]);
        printf("Ingrese el valor inicial de x1: ");
        scanf("%lf", &x0_newton[1]);
        double tolerancia_newton;
        printf("Ingrese la tolerancia: ");
        scanf("%lf", &tolerancia_newton);
        double h_newton;
        printf("Ingrese el valor de h: ");
        scanf("%lf", &h_newton);
        // creamos las matrices
        double **J2x2 = generar_matriz(2);
        double *F2 = generar_vector(2);
        // calculamos el método de newton
        newton_nolineal(2, J2x2, F2, x0_newton, tolerancia_newton, h_newton, N_newton, f1, f2);
        printf("El valor de x1 es: %lf\n", x0_newton[0]);
        printf("El valor de x2 es: %lf\n", x0_newton[1]);
        liberar_matriz(J2x2, 2);
        liberar_vector(F2);
        break;
    case 23:
        // Calcular Método de Broyden
        double x_broyden[2];
        printf("Ingrese el valor inicial de x0: ");
        scanf("%lf", &x_broyden[0]);
        printf("Ingrese el valor inicial de x1: ");
        scanf("%lf", &x_broyden[1]);
        double tolerancia_broyden;
        printf("Ingrese la tolerancia: ");
        scanf("%lf", &tolerancia_broyden);
        int iteraciones_broyden;
        printf("Ingrese la cantidad de iteraciones: ");
        scanf("%d", &iteraciones_broyden);
        broyden(x_broyden, tolerancia_broyden, iteraciones_broyden, f1, f2);
        printf("El valor de x1 es: %lf\n", x_broyden[0]);
        printf("El valor de x2 es: %lf\n", x_broyden[1]);
        break;
    case 24:
        // Calcular Gradiente Conjugado No Lineal
        int N_grad_conj;
        printf("Ingrese la cantidad de iteraciones: ");
        scanf("%d", &N_grad_conj);
        double x0_grad_conj[2];
        printf("Ingrese el valor inicial de x0: ");
        scanf("%lf", &x0_grad_conj[0]);
        printf("Ingrese el valor inicial de x1: ");
        scanf("%lf", &x0_grad_conj[1]);
        double tolerancia_grad_conj;
        printf("Ingrese la tolerancia: ");
        scanf("%lf", &tolerancia_grad_conj);
        gradiente_conjugado_nolineal(x0_grad_conj, tolerancia_grad_conj, N_grad_conj, f1, f2);
        printf("El valor de x1 es: %lf\n", x0_grad_conj[0]);
        printf("El valor de x2 es: %lf\n", x0_grad_conj[1]);
        break;
    case 25:
        printf("Saliendo del programa.\n");
        break;
    default:
        printf("Opción no válida.\n");
        break;
    }

    return 0;
}
