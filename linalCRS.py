# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 16:14:18 2021

@author: jor_a
"""

#Biblioteca de funciones para solucionar sistemas de ecuaciones lineales
#Y utilizar CRS

import numpy as np

def mat2crs(matriz):
    ###CRS##
    #Al utilizar python se considera el indice inicial = 0
    
    x,y = matriz.shape#se obtienen las dimensiones de la matriz
    val = []#lista para almacenar los valores no nulos
    col_ind = []#lista para almacenar los índices de la columna de cada valor
    ren_elem = []#Lista para almacenar la informacion de los renglones de la matriz
    cont = 0#indice inicial del vector de renglones
    
    for i in range(x):
        
        ren_elem.append(cont)
        
        for j in range(y):
            
                if matriz[i,j] != 0:
                    
                    val.append(matriz[i,j])#se alnacenan los valores de la matriz
                    col_ind.append(j)#se almacenan los valores de la columna
                    cont = cont+1#incremento de los indices de renglones
   
    #se transforman las listas a arreglos de numpy
    valor = np.array(val)
    col = np.array(col_ind)
    ren = np.array(ren_elem)
    
    return(valor,col,ren)

def maxerror(r):
    val = abs(r[0])
    for i in r:
        if abs(i) >= val:
            val = abs(i)

    return val


def prodmatCRS(valor,col,ren,vector):##Producto matricial utilizando CRS
    
    res=np.zeros(ren.size)
    
    for i in range(ren.size):
        suma = 0
        if i != ren.size-1:
            for j in range(ren[i],ren[i+1]):
                
                suma += valor[j]*vector[col[j]]
            
            res[i] = suma
        else:
            for j in range(ren[i],valor.size):
                suma += valor[j]*vector[col[j]]               
            res[i] = suma
    return res 
    
def prodpunto(vec1,vec2):
    if vec1.shape == vec2.shape:#los vectores deben tener el mismo numero de elementos
        suma = 0
        for i in range(vec1.shape[0]):
            
            suma += vec1[i]*vec2[i]
    return suma

def gradbic(mat,mat_trans,vector):
    
    dim = vector.shape[0]#tamaño de la matriz
    tol = 1e-5
    
    #Datos en CRS
    val,col_ind,ren_in =  mat[0], mat[1], mat[2]
    valt,col_indt,ren_int = mat[0], mat[1], mat[2]#Matriz transpuesta convertida a CRS
    
    
    #vectores iniciales
    x = np.ones(dim)#Solucion inicial, vector lleno de valor = 1     
    r = np.copy(vector)-prodmatCRS(val,col_ind,ren_in,x)#residuo(Gradiente inicial)
    r_ps = np.copy(r)# Pseudo gradiente inicial
    p = np.copy(r)#dirección de descenso inicial
    p_ps = np.copy(r_ps)
    
    maxitera = 1000
    i = 0
   
    while i < maxitera:
        if i == maxitera-5:
            print('no solucion')
        beta = prodpunto(r,r_ps)
        
        w = prodmatCRS(val,col_ind,ren_in,p)
        
        w_ps = prodmatCRS(valt,col_indt,ren_int,p_ps)
        
        
        rho = beta/(prodpunto(w, p_ps))
        
        x += rho * p
        r -= rho * w
        r_ps -= rho * w_ps
        #print('r', r)
        #print('r_ps', r_ps)
        
                
        gamma = beta
        beta = prodpunto(r_ps,r)
        
        alfa = (beta/gamma)
        
        p = r + alfa * p
        p_ps = r_ps + alfa * p_ps
        
        
        #comprobacion
        
        if maxerror(r)<=tol:#norma del error maximo
            
            #print(tol)
            #print(maxerror(r))
            break
        i += 1
    return(x)

   