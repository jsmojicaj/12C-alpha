## Codigos

Esta carpeta contiene los codigos creados que requieren de los datos experimentales, archivos txt guardados en la carpeta de [Datos](https://github.com/krishnamoji/12C-alpha/tree/main/Codigos/Datos), y genera una imagen que se guarda en la carpeta [Imagenes](https://github.com/krishnamoji/12C-alpha/tree/main/Codigos/Imagenes).

Se puede ejecutar desde el notebook [SeccionNuclear](https://github.com/krishnamoji/12C-alpha/blob/main/Codigos/SeccionNuclear_alfa_12C.ipynb) o ejecutando el archivo python [seccionnuclear](https://github.com/krishnamoji/12C-alpha/blob/main/Codigos/seccionnuclear.py). Se recomienda ejecutar el notebook por ser mas interactivo.

Para ejecutar el codigo primero se deben cargar los datos experimentales y después se debe ejecutar la funcion **Function** 

```python
Datos12C = np.loadtxt(DireccionDatos)[:,0:2]
Function(Ab=12, Zb=6, Ap=4, Zp=2, E=139, parametros=[108.1, 1.22, 0.76, 16.9, 1.85, 0.47, 1.26 ], Datos=Datos12C, ylim=(0.01,100), element='C')
```
Donde **A** y **Z** son los números de masa y atomicos, el indice **n** y **p** hacen referencia a si son del balnco o del proyectil; **E** es la energía de colisión; **parametros** es la lista de los coeficientes de ajuste del modelo óptico (potencial complejo de tipo Wood-Saxon); **Datos** son los datos experimentales extraidos de [NRV](http://nrv.jinr.ru/nrv/webnrv/expdata/?tab=elastic); **ylim** es el rango en el que se grafica la función; **element** es el simbolo del elemento objetivo.

<h3> Ejecución Local </h3>

Primero se descarga el repositorio desde la página del repositorio, [archivo comprimido](https://github.com/krishnamoji/12C-alpha/archive/refs/heads/main.zip), o desde la linea de comando de linux:

```linux
git init
git clone https://github.com/krishnamoji/12C-alpha.git
```
Para ejecutar el código ejecute el siguiente comando en una temrinal de linux

```linux
cd 12C-alpha/Codigos/
pip install -r librerias.txt
python3 SeccionNuclear_alfa_12C.py
```

Al ejecutar estos codigos se generarán varias imagnes en la carpeta de Imagenes contenida en Codigos. Otra forma de ejecutar el comando es mediante el archivo de JupyterNoteBook [SeccionNuclear_alfa_12C](https://github.com/krishnamoji/12C-alpha/blob/main/Codigos/SeccionNuclear_alfa_12C.ipynb) donde solo es dar click en ejecutar todas las celdas.

