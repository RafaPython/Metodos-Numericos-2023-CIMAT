
# Tarea 8: Método de Raleigh, QR y Subespacio
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en Linux (Fedora 38)

Menú de códigos para la resolución mediante el método de Rayleigh (eigenvector y eigenvalor), factorización QR, Método de iteración subespacio y método de Gradiente
conjugado y gradiente conjugado precondicionado para resolver sistemas de ecuaciones lineales. 

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
IMPORTANTE: Para el caso de los métodos de Gradiente Conjugado y Gradiente Conjugado Precondicionado, se pedirá hasta después de haberse seleccionado la opción 
el nombre del archivo del vector.

## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando dos métodos

1. Método de Iteración Subespacio
2. Método de Rayleigh
3. Método de factorización QR
4. Método de Gradiente Conjugado 
5. Método de Gradiente Conjugado Precondicionado

En cada método se verifica el tamaño de la matriz a abrir, si es menor a 10 se imprime en la consola, si es mayor a 10 se guarda en un archivo, se le pedirá al usuario ingresar el nombre del archivo con la extención que desee. 

### Uso y entrada

Ambos programas requieren del uso de la paquetería **_solvers.h_** dentro de la misma carpeta donde se encuentre el programa. Para utilizar el programa principal, debes proporcionar los siguientes argumentos (banderas) de línea de comandos:

- `NombreMatriz.txt`: Nombre del archivo que contiene la matriz del sistema.
- `TOLERANCIA`: (**_Double_**) Valor de tolerancia que se utilizará a lo largo de todos los cálculos del programa.
- `NUMERO_ITERACIONES`: (**_int_**) El número máximo de iteraciones (**IMPORTANTE**: entre mayor sea el tamaño de la matriz naturalmente más iteraciones tomará, por lo que se tendrá que tomar esto en cuenta al momento de ejecutar el programa)

El programa te dará a escoger cuál de los dos métodos utilizar. Si se elige el método del gradiente conjugado o el método del gradiente conjugado precondicionado este pedirá además 

- `Nombre_archivo_vector` : Nombre del archivo que contiene el vector para el sistema.

Después, el programa realizará el método que se eligió. 

### Salida
Para un archivo con dimensiones menores a 10, el programa imprime las matrices, los valores y vectores propios en la consola, para el caso de la factorización 
QR se imprime la matriz original y las matrices Q y R, asimismo se verifica que A - QR = 0. Para el caso del gradiente conjugado y gradietne conjugado precondicionado, se imprimen las soluciones al sistema. Para archivos de dimensiones mayores a 10, todo lo que anteriormente se imprimía en pantalla (a exepción de la matriz original y el vector original según sea el caso) se guardarán dentro de archivos .txt para lo cuál el usuario tendrá que administrar el nombre con el que desea guardarlos.
Los archivos generados se encontrarán dentro de la misma carperta donde se encuentre el programa.

NOTA: por alguna razón cuando estaba debuggeando con valgrind, comencé (por que al inicio no tenía) a tener problemas al momento de imprimir 
mis matrices, solamente con ello, desconozco si haya sido algún error mío o en mi computadora, no me presenta fugas de memoria pero me dice 
que estoy imprimiendo una matriz no inicializada por alguna razón, aunque pues claramente se puede imprimir sin problemas. Este problema se 
presenta en particular al momento de imprimir la matriz R en la factorización QR, pero la matriz Q se imprime sin problemas. 
