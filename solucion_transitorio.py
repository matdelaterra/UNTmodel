# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 13:40:51 2020

@author: jor_a
"""
#Solución para la ecuación de transporte con flujo en estado estacionario
#Utilizando diferencias finitas centradas en el tiempo, con el método de 
#diferencias hacia atrás en el tiempo usando upwind

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


q = 0.24#m/d
porosidad = 0.25
vz = q#/porosidad #cm/d
R = 1 #Retardación
dt = 1#dia
dz = 10 #cm
prof = 1000
L = 100
T = 2000#pasos
D = 0.001 #dispersión
dominio = np.linspace(0, L*dz, L+1)
dom = np.linspace(0, prof, int(prof//dz)+1)
#upwind
Dh = D*(1 + (vz*dz)/(2*D))

##Fuente/reacción
reac = 0#-0.1
#Condiciones iniciales
condicion = np.zeros(L+1)
frontera = 60
condicion[0] = frontera

incognitas = L


a = (Dh * dt)/(R*dz**2)
b = (vz * dt)/(R*2 * dz)

## tamaño del dominio
#pasos de tiempo

soluciones = [np.array([condicion])]

for t in range(T):#ciclo de tiempo
    matriz = np.zeros((incognitas, incognitas))
    vector = np.zeros(incognitas)
    for x in range(L):
        if x == 0:
            matriz[x,x] = 1 + 2*a
            matriz[x,x+1] = -a + b
                
            vector[x] =  condicion[x+1] + (a+b)*frontera + reac
                
        elif 0 < x < L-1:
            matriz[x,x+1] = -a + b
            matriz[x,x] = 1 + 2*a
            matriz[x,x-1] = -(a+b)
                
            vector[x] = condicion[x+1] + reac
        
        elif x == L-1:
            matriz[x,x-1] = -(a+b)
            matriz[x,x] = 1+a+b
                
            vector[x] = condicion[x+1] + reac
    
    sol = np.linalg.solve(matriz, vector)
    condicion[1:] = sol[:]
    soluciones.append(np.array([condicion]))


plt.plot(dominio, soluciones[-1].T, label='Upwind')


# fps = 30
# nSeconds = 5
# snapshots = [ np.random.rand(1,5) for _ in range( nSeconds * fps ) ]

# # First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()
ax.set(ylabel = 'Dominio Z', title='Concentración') 

im = ax.imshow(soluciones[0].T, aspect='auto', cmap='rainbow')
ax.xaxis.set_ticklabels('')
cbar = fig.colorbar(im, spacing='proportional',
                    shrink=0.9, ax=ax)
cbar.set_label("Concentración $mg/L$")

#plt.axis('off')
def animate_func(i):
    
    im.set_array(soluciones[i].T)
    ax.set
    return [im]

anim = animation.FuncAnimation(fig, 
                               animate_func,#, 
                               frames = len(soluciones),
                               repeat=False
                                #interval = 1000 / fps, # in ms
                                )

