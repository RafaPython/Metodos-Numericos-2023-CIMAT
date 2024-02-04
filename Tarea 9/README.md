
# Tarea 9: Métodos de interpolación
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en Linux (Fedora 38)

Menú de codigos para la resolución de valores interpolados mediante Taylor, Lagrange, Neville y Newton.

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
Para esta tarea no se creo un sistema con banderas para lectura de datos, solamente se le pide al usuario introducir toda la información necesaria
mediante la consola.

## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando dos métodos

1. Método de interpolación de Taylor
2. Método de interpolación de Lagrange
3. Método de interpolación de Neville
4. Método de interpolación de Newton

Cada método pide al usario introducir la información necesario así como los datos con los que va a interpolar.

### Uso y entrada
Cada método pide al usario introducir la información necesario así como los datos con los que va a interpolar.

### Salida
Se imprimirá el valor interpolado así como todos los datos ingresados dentro de la consola. 
