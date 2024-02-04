
# Tarea 12: Diferenciaci ́on y minimizaci ́on de funciones
### Nombre: Rafael Alejandro García Ramírez 
### email: rafael.ramirez@cimat.mx

Hecho en macOS (Ventura 13.3)

Menú de códigos para la resolución de diferenciación de funciones y datos y métodos para resolución de sistemas de ecuaciones no lineales por los métodos de punto fijo, Newton, Broyden y gradiente conjugado de Fletcher-Reeves.

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
NOTA: Para cada caso se pedirá introducir los valores de de interés del sistema, esto es en el caso de las derivadas de funciones (que se realizará para una función genérica que se puede moficar) serán los datos a evaluar y $h$. Para el caso de las derivadas de datos se pedirá los valores $f(x + ih)$ con $i \in \mathbb{N}_{0}$ del sistema. Asi mismo, para la resulución de un sistema de ecuaciones no lineales se tendrá un ejemplo genérico (el del problema 7a) con dos variables, las cuales también son modificables. Para el caso del jacobiano y hessiana se tienen 3 ecuaciones genéricas las cuales alimentan cada escenario (dos y tres variables).
## Descripción

El programa principal permite la resolución de sistemas de ecuaciones lineales utilizando dos métodos

1. Calcular derivada hacia adelante
2. Calcular derivada hacia atrás
3. Calcular derivada centrada
4. Calcular derivada de tres puntos en el extremo
5. Calcular derivada de tres puntos en el punto medio
6. Calcular derivada de cinco puntos en el punto medio
7. Calcular derivada de cinco puntos en el extremo
8. Calcular segunda derivada
9. Calcular derivada hacia adelante con datos
10. Calcular derivada hacia atrás con datos
11. Calcular derivada centrada con datos
12. Calcular derivada de tres puntos en el extremo con datos
13. Calcular derivada de tres puntos en el punto medio con datos
14. Calcular derivada de cinco puntos en el punto medio con datos
15. Calcular derivada de cinco puntos en el extremo con datos
16. Calcular segunda derivada con datos
17. Calcular Jacobiano 3x3
18. Calcular Jacobiano 2x2
19. Calcular Hessiana 3x3
20. Calcular Hessiana 2x2
21. Calcular Punto Fijo No Lineal
22. Calcular Método de Newton No Lineal
23. Calcular Método de Broyden
24. Calcular Gradiente Conjugado No Lineal

Cada método pide al usuario introducir la información necesario y con funciones genéricas en cada caso necesario.

### Uso y entrada
Cada método pide al usuario introducir la información/datos necesario para el cálculo sobre datos o funciones genericas segúns sea el caso.

### Salida
Se imprimirá la matriz/soluciones/derivada en pantalla según sea el caso.
