## Codigos

Esta carpeta contiene los codigos creados que requieren de los datos experimentales, archivos txt guardados en la carpeta de [Datos](https://github.com/krishnamoji/12C-alpha/tree/main/Codigos/Datos), y genera una imagen que se guarda en la carpeta [Imagenes](https://github.com/krishnamoji/12C-alpha/tree/main/Codigos/Imagenes).

Se puede ejecutar desde el notebook [SeccionNuclear](https://github.com/krishnamoji/12C-alpha/blob/main/Codigos/SeccionNuclear_alfa_12C.ipynb) o ejecutando el archivo python [seccionnuclear](https://github.com/krishnamoji/12C-alpha/blob/main/Codigos/seccionnuclear.py). Se recomienda ejecutar el notebook pues este es más interactivo y permite correr varios sistemas al tiempo.

Para ejecutar el codigo primero se deben cargar los datos experimentales y después se debe ejecutar la funcion **Function** 

[Descargar repositorio](https://github.com/krishnamoji/12C-alpha/archive/refs/heads/main.zip)

```python
Datos12C = np.loadtxt(DireccionDatos)[:,0:2]
Function(Ab=12, Zb=6, Ap=4, Zp=2, E=139, parametros=[108.1, 1.22, 0.76, 16.9, 1.85, 0.47, 1.26 ], Datos=Datos12C, ylim=(0.01,100), element='C')
```
Donde **A** y **Z** son los números de masa y atomicos, el indice **n** y **p** hacen referencia a si son del balnco o del proyectil; **E** es la energía de colisión; **parametros** es la lista de los coeficientes de ajuste del modelo óptico (potencial complejo de tipo Wood-Saxon); **Datos** son los datos experimentales extraidos de [NRV](http://nrv.jinr.ru/nrv/webnrv/expdata/?tab=elastic); **ylim** es el rango en el que se grafica la función; **element** es el simbolo del elemento objetivo.
