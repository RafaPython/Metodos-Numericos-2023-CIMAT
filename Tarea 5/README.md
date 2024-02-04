NOTA PRIMER PROBLEMA 

Dentro del primer problema no pude liberar correctamente la memoria para la variable "image" que es la imagen de lectura con la 
función que se nos pasó. Cuando lo intento me da error Segmentation Fault, aunque según mis pruebas con Valgrind, se sigue liberando la memoria.
Es decir, aunque tenga el error Segmentation Fault, igual se libera completamente la memoria; desconozco por qué sea esto, igual dejé la parte que libera 
esta memoria como comentario por si se desea ignorar ese error. De igual forma se me dice en valgrind que se intenta liberar memoria a que ya fue liberada.

Además de esto, a mi me da un warning al querer utilizar la función gets en lugar de fgets, aunque el Dr. Hasimoto nos ha recomendado utilizar mejor gets. 
