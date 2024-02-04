
# Tarea 11: Mínimos Cuadrados e Integración
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en Linux (Fedora 38)

Menú de codigos para la resolución de mínimos cuadrados regulados (con kernels Polinomiales, trigonométricas y exponenciales) y los métodos de integración de 
Newton-Cotes y Cuadratura Gaussiana para las funciones $\sin{x}$, $x^{2}\ln{x}$ y $x^{2}e^{-x}$.
## Compilación

Para la compilación del archivo principal ```main.c``` se hace el comando 

```bash
$ gcc main.c -o main -lm
```

## Ejecución

Para ejecutar el programa principal, utiliza el siguiente comando con SIN banderas:

```bash
$ ./main 
```
NOTA: En el caso de utilizar el método de mínimos cuadrados regulados, el programa pedirá el nombre de los archivos de datos y los respectivas decisiones y datos 
del usuario, después pedirá ingresar el nombre con el que se desea guardar los datos.

## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando dos métodos

1. Método de Mínimos cuadrados regulados 
2. Método de Newton-Cotes
3. Método de Cuadratura Gaussiana

Cada método pide al usuario introducir la información necesario así como los nombres de los archivos de datos con los que va a realizar los mínimos cuadrados
regulados.

### Uso y entrada
Cada método pide al usuario introducir la información necesaria (como nombre de archivos a utilizar en el caso del método de mínimos cuadrados regulados) así como los datos con los que va a interpolar/calcular según sea necesario.

### Salida
Se imprimirá el valor interpolado/calculado de ser el caso, de otra forma guardará un archivo .txt con los datos de los pesos del proceso de mínimos cuadrados regulados.
