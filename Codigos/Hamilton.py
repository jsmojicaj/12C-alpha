import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg as lg

#------------------FUNCIONES
def Laguerre(n, x):
    if (n == 0):
        f = 1; d = 0
    else:
        f = 1 - x; fm1 = 1
    for i in range(2,n+1):
        fm2 = fm1; fm1 = f
        f = ((2*i-1-x)*fm1 - (i-1)*fm2)/i
    return f
def derivada(n, x):
    dx=10**-5
    return (Laguerre(n,x+dx)-Laguerre(n,x-dx))/(2*dx)
def xGaussLag(a, n):
    eps=np.exp(-15) # relative precision of zeros
    X = [0]*n
    W = [0]*n
    for i in range(0,n):
        if (i == 0): # initial approximation for zeros (Stroud & Secrest)
             xi = 3/(1+2.4*n) # 1st zero
        elif (i == 1):
            xi = 15/(1+2.5*n) + X[0] # 2nd zero
        else:
            f = (1/(i+1)+2.55)/1.9
            xi = (1+f)*X[i-1] - f*X[i-2] # recurrence
        f = 9*np.exp(99)
        while (abs(f) > eps*abs(xi)): # Newton-Raphson refinement
            f = Laguerre(n,xi)/derivada(n, xi)
            xi -= f
        X[i] = xi
        W[i] = np.exp(xi)/(xi*derivada(n, xi)*derivada(n, xi))
    for i in range(0,n): X[i] += a # scaling to interval [a,+inf)
    return (X,W)
def funcion(a,n,k,x):
    [X,W]=xGaussLag(a,n)
    return (-1)**(k+1)*np.sqrt(X[k])*(x/X[k])*(Laguerre(n,x)/(x-X[k]))*np.exp(-x/2)

#-----------------------CODIGO MAIN
a=0
n=50 #Grado del polinomio alto

C=1 # 1; h^2/2m = 1Ryd*a0**2
h=1.5
K=2 # 2; e^2 / 4*pi*epsilon=2*a0*Ryd
[X,W]=xGaussLag(a,n); #Se guarda los ceros en X[i] y los pesos en W[i]

#Crear las mtrices T, V y H
T=[]
for i in range (n):
    T.append([0]*n)
V=[]
for i in range (n):
    V.append([0]*n)    
H=[]
for i in range (n):
    H.append([0]*n)    

for i in range(0,n):
    for j in range(0,n):
        if(i==j):
            T[i][j]=(1/h**2)*C*((4+(4*n+2)*X[i]-X[i]**2)/(12*X[i]**2))
            V[i][j]=-K/(X[i]*h)
            H[i][j]=T[i][j]+V[i][j]
        else: #en el paper es (-1)**(i-j+1) que empiezan desde i=j=1 pero aqui pusimo por (-1)**(i-j) por que empezamos desde i=j=0
            T[i][j]=(1/h**2)*(-1)**(i-j)*C*(((X[i]+X[j]))/(np.sqrt(X[i]*X[j])*(X[i]-X[j])**2))
            V[i][j]=0
            H[i][j]=T[i][j]+V[i][j]


        
values, vectors = lg.eig(H) #Calculo de los autovalores y autovectores
#print(values)
ordenados=sorted(values) #ordenar los autovalores

print(ordenados) #mostrar los autovalores ordenados