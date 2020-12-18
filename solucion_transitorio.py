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


q = 10#cm/d
porosidad = 0.4
vz = q/porosidad #cm/d
#R = 0.5
dt = 0.01#dia
dz = 1 #cm
L = 50
T = 300#pasos
D = 0.01 #dispersión
dominio = np.linspace(0,L,L+1)

#upwind
Dh = D*(1 + (vz*dz)/(2*D))


#Condiciones iniciales
condicion = np.zeros(L+1)
condicion[0] = 60
frontera = 60
incognitas = L


a = ( Dh * dt)/(dz**2)
b = (vz * dt)/( 2 * dz)

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
                
            vector[x] =  condicion[x+1] + (a+b)*frontera 
                
        elif 0 < x < L-1:
            matriz[x,x+1] = -a + b
            matriz[x,x] = 1 + 2*a
            matriz[x,x-1] = -(a+b)
                
            vector[x] = condicion[x+1] 
        
        elif x == L-1:
            matriz[x,x-1] = -(a+b)
            matriz[x,x] = 1+a+b
                
            vector[x] = condicion[x+1]
    
    sol = np.linalg.solve(matriz, vector)
    condicion[1:] = sol[:]
    soluciones.append(np.array([condicion]))



# fps = 30
# nSeconds = 5
# snapshots = [ np.random.rand(1,5) for _ in range( nSeconds * fps ) ]

# # First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()
ax.set(xlabel = 'Dominio x', title='Concentración') 

a = soluciones[0]
im = ax.imshow(a, aspect='auto', cmap='rainbow')
ax.yaxis.set_ticklabels('')
cbar = fig.colorbar(im, spacing='proportional',
                    shrink=0.9, ax=ax)
cbar.set_label("Concentración $mg/L$")
plt.grid()
#plt.axis('off')
def animate_func(i):
    im.set_array(soluciones[i])
    return [im]

anim = animation.FuncAnimation(
                                fig, 
                                animate_func,#, 
                                frames = len(soluciones)
                                #interval = 1000 / fps, # in ms
                                )
