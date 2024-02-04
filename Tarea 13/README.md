
# Tarea 13: M étodos clásicos de solución de EDO’s
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en macOS (Ventura 13.3)

Menú de códigos para la resolución de la ecuación diferencial $y' = y$ y para la solución al sistema de ecuaciones diferenciales Lotka-Volterra de la forma 

$\left\{\begin{matrix}
x'(t) = 0.4x(t) - 0.018x(t)y(t)\\ 
y'(t) = -0.8y(t) +0.023x(t)y(t)
\end{matrix}\right.$

El menú incluye los métodos de Euler, Heun, Taylor de 2do orden y Runge - Kutta 4, todos para la resolución de ambos casos.
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
NOTA: Para cada caso se pedirá introducir los datos de condiciones iniciales y la cantidad de intervalos con los que se desea realizar la resolucíón de cada problema.

## Descripción

El programa principal permite la resolución de la ecuación diferencial $y' = y$ y del sistema de ecuaciones diferenciales Lotka-Volterra con las siguientes opciones

1. Metodo de Euler (1 ecuacion)
2. Metodo de Heun (1 ecuacion)
3. Metodo de Taylor de segundo orden (1 ecuacion)
4. Metodo de Runge - Kutta 4 (1 ecuacion)
5. Metodo de Euler (sistema de ecuaciones)
6. Metodo de Heun (sistema de ecuaciones)
7. Metodo de Taylor de segundo orden (sistema de ecuaciones)
8. Metodo de Runge - Kutta 4 (sistema de ecuaciones)

Cada método pide al usuario introducir la información necesaria (como condiciones inciiales) y el rango y número de intervalos de interés.

### Uso y entrada
Cada método pide al usuario introducir la información/datos necesario para la resolución.

### Salida
Se guardará dentro de un archivo $\texttt{.txt}$ los datos generados por el métodos, con cada paso $x +ih$ que se haya pedido.
