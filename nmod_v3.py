# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:31:28 2020

@author: Jorge Antonio Matías López
"""

## Modelo transporte de nitratos V3
'''
En esta versión se van a manejar los atributos del objeto dentro de un dataFrame de pandas
si es conveniente se comenzará con el desarrollo de una interfaz gráfica en Tkinter, así como 
la lógica para utilizar capas diferentes
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

##Orientado a objetos, el objeto es modelo 
class Modelo:
    def __init__(self, D=100, dz=1, propiedades=None, NH4=60, NO3=0.1):
        self.__NH4 = NH4
        self.__NO3 = NO3
        self.__profundidad = D 
        self.__incremento = dz   
        self.__propiedades = propiedades
        self.__vz = propiedades['HLR']/propiedades['qs'] # CTE
        self.__Td = self.__incremento/self.__vz #cte
        self.__df_salida = None
        self.__p_profundidad = None
        self.__p_presion =  None
        self.__p_humedad = None
        self.__cont_carbon = None
        self.__WFP = None
        self.__R = None
        self.__krmax = None
        self.__Vmax = None
        self.__mass_fl = None
        self.__fsw_nt = None
        self.__fsw_dnt = None
        self.__C_NH4 = None
        self.__C_NO3 = None
        self.__Ctotal = None
        self.__allow_graph = False
 
    
    def __perfil_z(self):### Perfil de profundidad
        profundidad = self.__profundidad
        incremento = self.__incremento
        prof = np.linspace(0, profundidad, int(profundidad//incremento)+1)
        self.__p_profundidad = prof
        self.__df_salida = pd.DataFrame(data = prof, columns = ['z'])
        return prof
    
    
    def __presion(self, Yo=0): ## Perfil de presión
        dz = self.__incremento
        psi = np.zeros(len(self.__p_profundidad))
        propiedades = self.__propiedades
        aG = propiedades['aG']
        #Yo = propiedades['Yo']
        Ks = propiedades['Ks']
        HLR = propiedades['HLR']
        
        for ind, val in enumerate(psi):
            psi[ind] = Yo
            log = (HLR/Ks) * math.exp(-aG*Yo) * (math.exp(aG*dz) - 1) + 1
            Yo = Yo - dz + (1/aG) * math.log(log)
        self.__p_presion = psi
        return psi
    
    
    def __tetha(self):####### Perfil de humedad
        presion = self.__p_presion
        aVG = self.__propiedades['aVG']
        qs = self.__propiedades['qs']
        qr = self.__propiedades['qr']
        n  = self.__propiedades['n']
        m  = 1 - 1/n
        #print(aVG, qs, qr,n,m, sep='\n')
        tetha = np.zeros(presion.shape[0])
        
        for ind, val in enumerate(presion):
            if val > 0:
                tetha[ind] = qs   
            else:
                ss = qr + (qs-qr)/ math.pow((1 + math.pow(abs(aVG*val),n)),m)
                tetha[ind] = ss
        self.__p_humedad = tetha
        return tetha


    def __fWFP(self):### intentar con map fWFP(self, hum)
        
        qs = self.__propiedades['qs']
        swp = self.__propiedades['swp']
        hum = self.__p_humedad
        WFP = np.zeros(len(hum))
        for i, tetha in enumerate(hum):
            sw = tetha/qs
            if sw > swp:
                WFP[i] = sw
            else:
                WFP[i] = swp
        self.__WFP = WFP
        return WFP
    
    
    def __fcont_c(self):
        ac =  self.__propiedades['ac']
        cont_c = np.zeros(len(self.__p_profundidad))
        for i, z in enumerate(self.__p_profundidad):
            cont_c[i] = math.exp(-ac*z)
        self.__cont_carbon = cont_c
        self.__df_salida['fz'] = np.flip(cont_c)
        return cont_c
    
    
    def __fswnt(self):
        WFP = self.__WFP
        fs  = self.__propiedades['fs']
        fwp = self.__propiedades['fwp']
        swp = self.__propiedades['swp']
        e2  = self.__propiedades['e2']	 # e2
        e3  = self.__propiedades['e3']
        sl  = self.__propiedades['sl']
        sh  = self.__propiedades['sh']
        
        fsw_nt = np.zeros(len(WFP))
        
        for i, wfp in enumerate(WFP):
            if sh < wfp <= 1:
                fsw_nt[i] = fs + (1-fs)*math.pow((1-wfp)/(1-sh),e2)
            elif sl <= wfp <= sh:
                fsw_nt[i] = 1
            elif swp <= wfp <= sl:
                fsw_nt[i] = fwp + (1-fwp)*math.pow((wfp-swp)/(sl-swp),e3)
        self.__fsw_nt = fsw_nt
        self.__df_salida['fsw_nt'] = np.flip(fsw_nt)
        return fsw_nt
        
      
    def __fswdnt(self):
        WFP  = self.__WFP
        sdn  = self.__propiedades['sdn']
        ednt = self.__propiedades['ednt']
        
        fsw_dnt = np.zeros(len(WFP))
        
        for i, wfp in enumerate(WFP):
            if wfp < sdn:
                fsw_dnt[i] = 0
            else:
                fsw_dnt[i] = math.pow((wfp-sdn)/(1-sdn),ednt)
        self.__fsw_dnt = fsw_dnt
        self.__df_salida['fsw_dnt'] = np.flip(fsw_dnt)
        return fsw_dnt
    
 ##### Parámetros que se utilizan en el calculo de C   
    def __fR(self):
        kd = self.__propiedades['kd']
        tetha = self.__p_humedad
        rho =  self.__propiedades['rho']
        r = np.zeros(len(tetha))
        
        for i, h in enumerate(tetha):
            r[i] = 1 + (rho*kd/(h))
        self.__R = r
        self.__df_salida['R'] = np.flip(r)
        return r
    
    def __f_krmax(self):
        kr_max = self.__propiedades['kr_max']
        kr = self.__fsw_nt * kr_max
        self.__krmax = kr
        self.__df_salida['krmax'] = np.flip(kr)
        return kr
    
    
    def __fVmax(self):
        Vmax = self.__propiedades['Vmax']
        Vm = self.__fsw_dnt * Vmax
        self.__Vmax = Vm
        self.__df_salida['Vmax'] = np.flip(Vm)
        return Vm
    
    
###################    
    def __fconc(self, km, Co, mu_max, R=1, fz=1):
        td = self.__Td
        ## Newton rhapson
        ini = 0.000001
        b = -Co-km*math.log(Co) + R*mu_max*fz*td
        con = lambda c,b, km : c + km * math.log(c) + b
        der_con = lambda c, km: 1 + km/c
        
        it = 1
        while True:
            nv = ini - con(ini,b,km)/der_con(ini,km)
            if abs(ini - nv) < 0.0000001 or it > 50:
                conc = nv
                break
            else:
                ini = nv
                it += 1
        return conc
    

    def __Cnit(self):
        #parametros
        C_ini_NH4 = self.__NH4
        C_ini_NO3 = self.__NO3
        R_nt = self.__df_salida['R']
        krmax = self.__df_salida['krmax']
        vmax = self.__df_salida['Vmax']       
        
        n = np.shape(self.__p_profundidad)[0] 
        c_NH4 = np.zeros(n)
        c_NH4[0] = C_ini_NH4
        c_NO3 = np.zeros(n)
        c_NO3[0] = C_ini_NO3
        
        for i in range(1,n):
            if C_ini_NH4 > 0.001:
                c_NH4[i] = self.__fconc(self.__propiedades['Km_nit'], C_ini_NH4, krmax[i], R_nt[i])
            else:
                c_NH4[i] = 0
            
            C_ini_NH4 = c_NH4[i]
            
            C_ini_NO3 = c_NO3[i-1] + c_NH4[i-1] - c_NH4[i]
    
            c_NO3[i] = self.__fconc( self.__propiedades['Km_dnt'],C_ini_NO3, vmax[i])
            
        self.__C_NH4 = c_NH4
        self.__df_salida['cNH4'] = c_NH4
        self.__C_NO3 = c_NO3
        self.__df_salida['cNO3'] = c_NO3
        self.__Ctotal = c_NH4 + c_NO3    
        self.__df_salida['Total N'] = self.__Ctotal
        
        #Flujo de masa
        Ctotal = self.__Ctotal
        hlr = self.__propiedades['HLR']
        mass_fl = np.zeros(Ctotal.shape[0])
        for i in range(Ctotal.shape[0]):
            mass_fl[i] = (2/3) * Ctotal[i]* (2/3) * hlr / 100
        
        self.__mass_fl = mass_fl
            
            
    def run(self):
        #cálculo de profundidad
        self.__perfil_z()
        #### HIDRAULICA Y CONTENIDO DE HUMEDAD
        #Cálculo de presiones
        self.__presion()# = 
        #perfil de humedad
        self.__tetha()
        #WFP
        self.__fWFP()
        #FACTORES PARA LAS REACCIONES
        ## efecto del contenido de carbono
        self.__fcont_c()
        ## efecto de sw en la nit
        self.__fswnt()
        ## efecto de sw en la desn
        self.__fswdnt()
        #ELEMENTOS PARA CALCULAR LAS CONCENTRACIONES
        ### factor de retardación
        self.__fR()
        ### kr max
        self.__f_krmax()
        ##Vmax
        self.__fVmax()
        ## Concentraciones
        self.__Cnit()
        
        self.__allow_graph = True
            
    def graficar(self, mass = False):
        if self.__allow_graph:
            
            x = self.__p_profundidad
            y1 = self.__C_NH4
            y2 = self.__C_NO3
            y3 = self.__Ctotal
           
            
            
            fig, ax = plt.subplots()
            ax.plot(x, y1)
            ax.plot(x, y2)
            ax.plot(x, y3)
            ax.set_xlabel('Profundidad')
            ax.set_ylabel('Concentración')
            ax.set_title('Concentración N')
            ax.legend(['CNH4','CNO3','CTotal'])
            
            if mass:
                y_fl = self.__mass_fl
                fig1, ax1 = plt.subplots()
                ax1.plot(x, y_fl)
                ax1.set_xlabel('Profundidad')
                ax1.set_ylabel('Mass flux $gm^2 d^{-1}$')
                ax1.set_title('Mass Flux TN')
                ax1.legend(['Mass Flux TN'])
                plt.show()
        else:
            print('Se debe ejecutar el modelo')
        


if __name__ == '__main__':
    #Diccionario de propiedades
    dic_h = {
    'HLR':5.00,
    'aG':0.025,#	 aG
    'aVG':0.015,#	aVG
    'Ks':14.75,#	Ks
    #'Yo':0,#	Yo
    'qr':0.0980,#	qr
    'qs':0.459,#	qs 
    'n':1.26,#	n 
    'l':0.50,	#  l 
    ######## nit
    'swp':0.15,# swp ### sirve para calcular WFP
    'fs':0.00,#	  fs 
    'fwp':0.00,#	  fwp
    'e2':2.27,	 # e2
    'e3':1.10,	#  e3 
    'kr_max':50.00,#	  kr max
    'Km_nit':5.00,#3	  Km,nit
    'bnit':0.35,#	  bnit 
    'sl':0.67,#	 sl
    'sh':0.81,#	  sh
    ####### dnt
    'ednt':3.77,#	  ednt
    'Vmax':2.56,#	  Vmax 
    'Km_dnt':5.00,#	  Km,dnt
    'bdnt':0.35,#	  bdnt 
    'sdn':0.00,#	  sdn 
    'ac':0.00,#	  ac 
    ##### R
    'kd':1.46,
    'rho':1.50	  #rho
    }

    #Creación del objeto modelo
    mod = Modelo(120, 1, propiedades=dic_h, NH4=10, NO3 = 1)
    mod.run()
    mod.graficar()


