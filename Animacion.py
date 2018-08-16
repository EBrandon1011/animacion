# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 08:46:12 2018

@author: brandon
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 08:37:32 2018

@author: brandon
"""
g=9.81  #Aceleración gravitacional.
k=32.0   #Constante de resorte [N/m].
br=0.0  #Constante de amortiguamiento del resorte.  br=$\beta$
bp=0.0  #Constante de amorgiguamiento del péndulo.  bp=$\gamma$
M=0.200   #Masa del objeto - Resorte.
m=0.96   #Masa del objeto - Péndulo.
l=0.5   #Longitud del péndulo.

from numpy import cos, sin, pi
import numpy as np
import matplotlib.pyplot as plt

#Definiendo las ecuaciones diferenciales dq_i/dt=f_i(q):

def dx(v):              #Dif. x --> dx/dt
    return v
    
def dv(x, v, ang, w):   #Dif. v --> dv/dt
    num=2*bp*cos(ang)*w-2*br*l*v+m*g*l*sin(ang)*cos(ang)-k*l*x+m*l**2*sin(ang)*w**2
    den=l*(M+m*(sin(ang))**2)
    return num/den
    
def dang(w):            #Dif. o --> do/dt
    return w
    
def dw(x, v, ang, w):   #Dif. w --> dw/dt
    num=2*bp*(m+M)*w-2*br*(m+M)*v+m*g*l*(m+M)*sin(ang)+m**2*l**2*sin(ang)*cos(ang)*w**2-m*l*k*cos(ang)*x
    den=m*l**2*(M+m*(sin(ang))**2) 
    return (-1)*num/den
    
#----------------------------------------------#   
#APLICACIÓN DEL MÉTODO DE EULER
    
def Euler(x0, v0, ang0, w0, T, NI):
    #Condiciones iniciales
    x=[]
    v=[]
    ang=[]
    w=[]
    time=[]
    x.append(x0)
    v.append(v0)   
    ang.append(ang0)
    w.append(w0)
    time.append(0)
    h=T/NI     #Tamaño de paso
    n=1
    while n<NI:
        pw=w[n-1]+h*dw(x[n-1],v[n-1],ang[n-1],w[n-1])   #Obteniendo  w[n]
        pv=v[n-1]+h*dv(x[n-1],v[n-1],ang[n-1],w[n-1])   #Obteniendo  v[n]
        w.append(pw)
        v.append(pv)
        ang.append(ang[n-1]+h*dang(w[n-1]))
        x.append(x[n-1]+h*dx(v[n-1]))
        time.append(n*h)
        n+=1
    return x, v, ang, w, time

#Condiciones iniciales
x0=0.0
v0=0.0
ang0=170*pi/180
w0=0.0
T=20
N=T*1000
x,v,ang,w,t=Euler(x0, v0, ang0, w0, T, N)

#----------------------------------------------#

## POSICIÓN DEL PÉNDULO
xp=[]   #Posición del péndulo
yp=[]
yr=[]
for i in range(N):
    ypi=(-1)*l*cos(ang[i])
    xpi=x[i]+l*sin(ang[i])
    yr.append(0)
    xp.append(xpi)
    yp.append(ypi)


#----------------------------------------------#

## ANIMACIÓN DEL RESORTE Y PÉNDULO

i=0
while i<N:
    plt.xticks(np.arange(-1.0,1.25, step=0.25))
    plt.xlim(-1.1,1.1)
    plt.ylim(-0.7,0.7)
    plt.plot([-1.1,1.1],[0,0],'#A9A9A9')
    plt.plot([0,0],[-1.1,1.1],'#A9A9A9')
    plt.plot([-1.1,x[i]],[0,0],'#2F4F4F',linestyle='-.', linewidth=8)    #Resorte
    plt.plot(x[i],yr[i],'#191970',marker='s', markersize=12)             #Carrito
    plt.plot(xp[i],yp[i],'#191970', marker='o', markersize=12)           #Péndulo
    x1,y1=x[i],yr[i]
    x2,y2=xp[i],yp[i]
    plt.plot([x1,x2],[y1,y2],'#2F4F4F',linewidth=2)                     #Hilo del péndulo
    plt.savefig('170gPR{}.png'.format(i),dpi=200)
    #plt.show()
    plt.cla()
    i+=10
