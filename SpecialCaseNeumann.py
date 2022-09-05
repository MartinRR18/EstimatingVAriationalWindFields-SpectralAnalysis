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
enesima = 10

#Define the function f(x) as executable 
def f(z):
    f = z*exp(-z/H)
    return f

def dg(x):
    dg = ((np.pi)/(2*a))*sin((np.pi)/(2*a)*x)
    return dg

def e(y):
    e = 1
    return e

def F(x,y,z):
    F = dg(x)*e(y)*f(z)
    return F


#Base con condiciones Neumann normalizada; componente x
varphi_i=np.sqrt(2/a)*cos((m*np.pi*x)/a)

#Base con condiciones Neumann normalizada; componente y
varphi_j=np.sqrt(2/b)*cos((n*np.pi*y)/b)

#Base con condiciones Neumann normalizada; componente z
varphi_k=np.sqrt(2/c)*cos((l*np.pi*z)/c)



#------------------Calculamos los coeficientes gi y ej----------------
#-----------considerando condiciones Neumann en las laterales-------
dgi1=integrate(dg(x)*varphi_i,(x,0,a))

ej1=integrate(e(y)*varphi_j,(y,0,b))


#-----------------Calculamos los coeficientes fk----------------------
#-----------considerando condiciones Mixtas en la vertical------------
#---........Establacemos parámetros necesarios--------
w1=(l*np.pi/c)
gama=H/(1+((w1*H)**2))

fk1=((gama)**2)*np.sqrt(2/c)*(((-1)**l)*exp(-c/H)*((c/gama)-1+((H*w1)**2))+1-((H*w1)**2))


#----------------------------SERIES&PLOTS------------------------------

#inf=float(integrate((1/a)*np.sqrt(2/a)*dg(x)**2,(x,0,a)))
#inic = np.sqrt(inf)

dg0=integrate(dg(x),(x,0,a))
serie_dgm=(1/a)*dg0
for i in range(1,50):
    serie_dgm = serie_dgm + ((varphi_i*dgi1).subs(m,i))    
plot(serie_dgm, xlim=(-2.5, 20.5), ylim=(-1.5,1.5)) 

e0 = integrate((1/b)*e(y),(y,0,b))
serie_en = e0
for j in range(1,30):
    serie_en = serie_en + ((varphi_j*ej1).subs(n,j))
#Usando el modulo para graficas de sympy
plot(serie_en, xlim=(-2.5, 10.5), ylim=(-1.5,1.5)) 

f0 = integrate((1/c)*f(z),(z,0,c))
serie_fl = f0
for k in range(1,30):
    serie_fl = serie_fl + ((varphi_k*fk1).subs(l,k))
plot(serie_fl, xlim=(-1.5, 5.5), ylim=(-1.5,1.5))#(z,-6,6))


#----------------------Norma de F(x,y,z)------------------------------

fc = integrate(F(x,y,z)**2,(x,0,a),(y,0,b),(z,0,c))
normaF = np.sqrt(float(fc))
print( 'La norma de F(x,y,z)=',F(x,y,z),' en la región xM=yM=20km por zM=5km es: ', normaF)
