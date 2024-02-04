
# Tarea 6: Método de Jacobi y Gauss-Seidel
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en Linux (Fedora 38)

Resuelve sistemas de ecuaciones lineales utilizando el método de Jacobi y el método De Gauss-Seidel.

## Compilación

Para la compilación del archivo principal ```main.c``` se hace el comando 

```bash
$ gcc main.c -o main -lm
```

## Ejecución

Para ejecutar el programa principal, utiliza el siguiente comando con las banderas:

```bash
$ ./main NombreMatriz.txt NombreVector.txt TOLERANCIA NUMERO_ITERACIONES
```


## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando dos métodos

1. Método de Jacobi
2. Método de Gauss-Seidel

En cada método se verifica el tamaño de la matriz a abrir, si es menor a 10 se imprime en la consola, si es mayor a 10 se guarda en un archivo con el nombre del método más un identificador del tamaño de la matriz (en formato .txt)

### Uso y entrada

Ambos programas requieren del uso de la paquetería **_solvers.h_** dentro de la misma carpeta donde se encuentre el programa. Para utilizar el programa principal, debes proporcionar los siguientes argumentos (banderas) de línea de comandos:

- `NombreMatriz.txt`: Nombre del archivo que contiene la matriz del sistema.
- `NombreVector.txt`: Nombre del archivo que contiene el vector de términos independientes.
- `TOLERANCIA`: (**_Double_**) Valor de tolerancia que se utilizará a lo largo de todos los cálculos del programa.
- `NUMERO_ITERACIONES`: (**_int_**) El número máximo de iteraciones (**IMPORTANTE**: entre mayor sea el tamaño de la matriz naturalmente más iteraciones tomará, por lo que se tendrá que tomar esto en cuenta al momento de ejecutar el programa)

El programa te dará a ecoger cuál de los dos métodos utilizar. Después, el programa realizará los siguientes pasos:

1. Verificación de la condición de diagonalmente dominante, en caso de no serlo te dará a escoger entre continuar con el método o abandonar el programa.
2. Inicialización del vector $x_{1}$ con números 1
3. Iteración Jacobi/Gauss-Seidel
4. Verificación de condición de convergencia
5. Actualización de valores de $x_{0}$ y $x_{1}$.
6. Verificación de $Ax-b = 0$

### Salida

Para un archivo de dimensiones menores a 10, el programa imprime $A$,$b$ y $x_{0}$ inicial para después imprimir la cantidad de iteraciones que se llevaron a cabo y la solución $x_{0}$ final, en caso contrario de tratarse de dimensiones mayores a 10, el programa solamente genera archivos con los resultados de la solución e indica la cantidad de iteraciones llevadas a cabo. 

En ambos casos se verificará $Ax-b= 0$ y se imprimirá la cantidad de elementos $x_{j}$ que no cumplieron con esta condición.

Los archivos generados se encontrarán dentro de la misma carperta donde se encuentre el programa.
