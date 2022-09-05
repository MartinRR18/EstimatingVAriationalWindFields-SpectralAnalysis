#Import all sympy module; enable us to work the equations 
#in the symbolic way. 
from sympy import *
#Import also the symbolic variables 'n' and 'x'
from sympy.abc import x, y, z, n, m, l, a
import numpy as np
import matplotlib.pyplot as plt


#Set a,b,c as parameters for limit of integration 
#and beta as mean velocity value with dimensions of [km*s^-1].

a=31.2
b=31.2
c=4
beta=0.00001
H=10
#Set the number of coeficients in the Fourier series.
#enesima = 10

#Define the function f(x) as executable 
#def f(z):
f = z*exp(-z*0.01)
#    return f

#def dg(x):
dg = ((np.pi)/(2*a))*sin((np.pi)/(2*a)*x)
#    return dg

#def e(y):
e = 1
#    return e

def F(x,y,z):
    F = dg*e*f
    return F


#Base de condiciones Dirichlet normalizada; componente x
varphi_i=np.sqrt(2/a)*sin((m*np.pi*x)/a)

#Base de condiciones Dirichlet normalizada; componente y
varphi_j=np.sqrt(2/b)*sin((n*np.pi*y)/b)

#Base de condiciones Mixtas normalizada; componente z
varphi_k=np.sqrt(2/c)*cos(((l+(1/2))*np.pi*z)/c)



#------------------Calculamos los coeficientes gi y ej----------------
#-----------considerando condiciones dirichlet en las laterales-------
dgi1=integrate(dg*varphi_i,(x,0,a))

ej1=integrate(e*varphi_j,(y,0,b))


#-----------------Calculamos los coeficientes fk----------------------
#-----------considerando condiciones Mixtas en la vertical------------
fk1=integrate(f*varphi_k,(z,0,c))


#-----------------Coeficiente fk----------------------
#------considerando condiciones Mixtas en la vertical------------
#---........Establacemos parámetros necesarios--------
H=10 #km
w2=((1/2)+l)*(np.pi/c)
gama=H/(1+(w2*H)**2)

#fk1=np.sqrt(2/c)*(((-1)**l)*gama*H*w2*exp(c/H)*(c+2*gama)-gama*(1+gama*H*w2))

#Coeficiente_corregido
fk1=np.sqrt(2/c)*gama*(((-1)**l)*c*H*w2*exp(-c/H)+gama*(1-(H*w2)**2))



#----------------------------SERIES&PLOTS------------------------------

serie_dgm = 0
for i in range(0,20):
    serie_dgm = serie_dgm + ((varphi_i*dgi1).subs(m,i))    
plot(serie_dgm,(x,-5,25))
           
serie_en = 0
for j in range(0,20):
    serie_en = serie_en + ((varphi_j*ej1).subs(n,j))
#Usando el modulo para graficas de sympy
plot(serie_en,(y,-5,25))

serie_fl = 0
for k in range(0,20):
    serie_fl = serie_fl + ((varphi_k*fk1).subs(l,k))
plot(serie_fl,(z,-6,6))



#-...-Límites de integración inicialización-------------------
xyM1 = 31.2
xym1 = 0
zM1 = 4
zm1 = 0

#-----------Calclulo de la integral de Fmnl2-----------------

int_ser1=integrate(serie_dgm,(x,xym1,xyM1))
int_ser2=integrate(serie_en,(y,xym1,xyM1))
int_ser3=integrate(serie_fl,(z,zm1,zM1))
int_Fmnl21 = int_ser1*int_ser2*int_ser3
    
#-------------Calclulo de la integral de F2------------------
int_F21 = integrate(F(x,y,z),(x,xym1,xyM1),(y,xym1,xyM1),(z,zm1,zM1))
flux_11 = beta*(int_F21-int_Fmnl21)
    
#------- ---Porcentaje de Masa que fluye en 3hr---------------

PorcentajeMasa1=(3*100*3600*abs(flux_11))/((xyM1-xym1)*(xyM1-xym1)*(zM1-zm1))
    
#---------------------Resultados------------------------------
print('\n\n Flujo en la subregión (',xym1,',',xyM1,')^{2}x(',zm1,',',zM1,') --> ',flux_11)
print(' %Masa que fluye en la subregión en 3hr   --> ',PorcentajeMasa1)


