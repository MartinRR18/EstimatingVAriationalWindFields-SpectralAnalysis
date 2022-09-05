#Import all sympy module; enable us to work the equations 
#in the symbolic way. 
from sympy import *
#Import also the symbolic variables 'n' and 'x'
from sympy.abc import x, y, z, n, m, l, a
import numpy as np
import matplotlib.pyplot as plt


#Set a,b,c as parameters for limit of integration 
#and beta as mean velocity value with dimensions of [km*s^-1].

a=20
b=20
c=5
beta=0.001
H=10 #km
#Set the number of coeficients in the Fourier series.
#enesima = 10

#Define the function f(x) as executable 
def f(z):
    f = z*exp(-z/H)
    return f

def dg(x):
    dg = -((np.pi)/(2*a))*sin(((np.pi)/(2*a))*x)
    return dg

def e(y):
    e = 1
    return e

def F(x,y,z):
    F = dg(x)*e(y)*f(z)
    return F


#Base de condiciones Dirichlet normalizada; componente x
varphi_i=np.sqrt(2/a)*sin((m*np.pi*x)/a)

#Base de condiciones Dirichlet normalizada; componente y
varphi_j=np.sqrt(2/b)*sin((n*np.pi*y)/b)

#Base de condiciones Mixtas normalizada; componente z
varphi_k=np.sqrt(2/c)*cos(((l+(1/2))*np.pi*z)/c)



#------------------Calculamos los coeficientes gi y ej----------------
#-----------considerando condiciones dirichlet en las laterales-------
dgi1=integrate(dg(x)*varphi_i,(x,0,a))

ej1=integrate(e(y)*varphi_j,(y,0,b))


#-----------------Coeficiente fk----------------------
#------considerando condiciones Mixtas en la vertical------------
#---........Establacemos parámetros necesarios--------
w2=((1/2)+l)*(np.pi/c)
gama=H/(1+(w2*H)**2)

fk1=np.sqrt(2/c)*gama*(((-1)**l)*c*H*w2*exp(-c/H)+gama*(1-H*w2))




#------------------------SERIES&PLOTS------------------------------

serie_dgm = 0
for i in range(0,30):
    serie_dgm = serie_dgm + ((varphi_i*dgi1).subs(m,i))    
plot(serie_dgm,(x,-5,25))
           
serie_en = 0
for j in range(0,30):
    serie_en = serie_en + ((varphi_j*ej1).subs(n,j))
#Usando el modulo para graficas de sympy
plot(serie_en,(y,-5,25))

serie_fl = 0
for k in range(0,30):
    serie_fl = serie_fl + ((varphi_k*fk1).subs(l,k))
plot(serie_fl,(z,-6,6))


#----------------------Norma de F(x,y,z)-------------------------

fc = integrate(F(x,y,z)**2,(x,0,a),(y,0,b),(z,0,c))
normaF = np.sqrt(float(fc))
print( 'La norma de F(x,y,z)=',F(x,y,z),' en la región xM=yM=20km por zM=5km es: ', normaF)

    
''' 
#----------------------Norma de Fmnl----------------------------
normaFmnl=0
for i in range(0,20):
    for j in range(0,20): 
        for k in range(0,40):
            normaFmnl = normaFmnl + (((dgi1)**2).subs(m,i))*(((ej1)**2).subs(n,j))*(((fk1)**2).subs(l,k))
    #En cada paso va sumando la norma término a termino        
print("La norma de la serie Fnml con m=n=l=20 es: ", sqrt(normaFmnl))
'''
