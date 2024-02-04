
# Tarea 7: Método Potencia y Potencia Inversa
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
$ ./main NombreMatriz.txt TOLERANCIA NUMERO_ITERACIONES
```
El número de iteraciones deberá ser sin notación científica. 

## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando dos métodos

1. Método de la Potencia
2. Método de la Potencia Inversa
3. Calcular los N valores propios más grandes/pequeños
4. Método de Jacobi para valores y vectores propios

En cada método se verifica el tamaño de la matriz a abrir, si es menor a 10 se imprime en la consola, si es mayor a 10 se guarda en un archivo con el nombre del método más un identificador del tamaño de la matriz (en formato .txt)

### Uso y entrada

Ambos programas requieren del uso de la paquetería **_solvers.h_** dentro de la misma carpeta donde se encuentre el programa. Para utilizar el programa principal, debes proporcionar los siguientes argumentos (banderas) de línea de comandos:

- `NombreMatriz.txt`: Nombre del archivo que contiene la matriz del sistema.
- `TOLERANCIA`: (**_Double_**) Valor de tolerancia que se utilizará a lo largo de todos los cálculos del programa.
- `NUMERO_ITERACIONES`: (**_int_**) El número máximo de iteraciones (**IMPORTANTE**: entre mayor sea el tamaño de la matriz naturalmente más iteraciones tomará, por lo que se tendrá que tomar esto en cuenta al momento de ejecutar el programa. NO introducirse en notación científica.)

El programa te dará a ecoger cuál de los dos métodos utilizar. Si se escoge el calcular N valores propios más grandes/pequeños, se pedirá especificar de qué tipo y después cuántos desea calcular.

### Salida

Para un archivo de dimensiones mayores o igual a 10, el programa guarda los valores propios y vectores propios dentro de un archivo el cuál tendrá en su primera fila datos respecto a la dimensión de la información contenida, esto es $\texttt{ancho}\times\texttt{alto}$. Para dimensiones menores a 10 el programa desplegará tanto los valores propios como los vectores propios en la consola.