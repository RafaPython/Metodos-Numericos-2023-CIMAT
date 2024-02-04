#include <stdio.h>
#include <math.h>

double e(double n){
    return pow((1 + (1/n)), n);
}

double e_mod(double n){
    return (1 / n);
}
int main() {
    const int k = 20;
    double datos[k - 1];
    
    for(int i = 0; i < k; i++) datos[i] = e(pow(10,i));
    
    printf("EL valor %d es = %.16f\n", k , datos[k - 1]);
    double errores[k - 1];
    double e1 =  exp(1);    
    
    for(int i = 0; i < k; i++) errores[i] = fabs(datos[i] - e1);
    printf("EL error en n = %d es: %.32f\n", k , errores[k - 1]);
    printf("\n");
    double peso[k - 1];
    for(int i = 0; i < k; i++) peso[i] = e_mod(pow(10,i));
    printf("El peso 1/n en n = %d es: %.32f\n", k , peso[k - 1]);
    // epsilon 
    
    printf("\nEPSILON\n");
    
    double epsilon = 1;
    int iterador = 0;
    while(1 + epsilon  > 1){
        epsilon /= 2;
        iterador++;
    }
    printf("El Epsilon en la iteracion %d es de: %.32f\n", iterador, epsilon);
    
    printf("\nLa relación (epsilon > 1/n) es: %s\n", (epsilon > peso[k - 1]) ? "True" : "False");
    printf("La diferencia entre epsilon y el ultimo 1/n es de: %.32f\n", fabs(epsilon - peso[k - 1]));
    printf("\n");

    return 0;
}