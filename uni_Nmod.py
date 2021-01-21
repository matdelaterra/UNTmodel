# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 18:25:45 2021

@author: jor_a
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

##Orientado a objetos, el objeto es modelo 
class Modelo:
    def __init__(self, D=100, dz=1, propiedades=None, NH4=60, NO3=0.1):
        self.NH4 = NH4
        self.NO3 = NO3
        self.profundidad = D 
        self.incremento = dz   
        self.propiedades = propiedades
        self.vz = propiedades['HLR']/propiedades['qs'] # CTE
        self.df_salida = None
        self.p_profundidad = None
        self.p_presion =  None
        self.p_humedad = None
        self.cont_carbon = None
        self.WFP = None
        self.R = None
        self.krmax = None
        self.Vmax = None
        self.fsw_nt = None
        self.fsw_dnt = None
        
    
    def perfil_z(self):### Perfil de profundidad
        profundidad = self.profundidad
        incremento = self.incremento
        prof = np.linspace(0, profundidad, int(profundidad//incremento)+1)
        self.p_profundidad = prof
        self.df_salida = pd.DataFrame(data = prof, columns = ['z'])
        return prof
    
    
    def calc_presion(self, Yo=0): ## Perfil de presión
        dz = self.incremento
        psi = np.zeros(len(self.p_profundidad))
        propiedades = self.propiedades
        aG = propiedades['aG']
        #Yo = propiedades['Yo']
        Ks = propiedades['Ks']
        HLR = propiedades['HLR']
        
        for ind, val in enumerate(psi):
            psi[ind] = Yo
            log = (HLR/Ks) * math.exp(-aG*Yo) * (math.exp(aG*dz) - 1) + 1
            Yo = Yo - dz + (1/aG) * math.log(log)
        self.p_presion = psi
        return psi
    
    
    def calc_humedad(self):####### Perfil de humedad
        presion = self.p_presion
        aVG = self.propiedades['aVG']
        qs = self.propiedades['qs']
        qr = self.propiedades['qr']
        n  = self.propiedades['n']
        m  = 1 - 1/n
        #print(aVG, qs, qr,n,m, sep='\n')
        tetha = np.zeros(presion.shape[0])
        
        for ind, val in enumerate(presion):
            if val > 0:
                tetha[ind] = qs   
            else:
                ss = qr + (qs-qr)/ math.pow((1 + math.pow(abs(aVG*val),n)),m)
                tetha[ind] = ss
        self.p_humedad = tetha
        return tetha


    def fWFP(self):### intentar con map fWFP(self, hum)
        
        qs = self.propiedades['qs']
        swp = self.propiedades['swp']
        hum = self.p_humedad
        WFP = np.zeros(len(hum))
        for i, tetha in enumerate(hum):
            sw = tetha/qs
            if sw > swp:
                WFP[i] = sw
            else:
                WFP[i] = swp
        self.WFP = WFP
        return WFP
    
    
    def fcont_c(self):
        ac =  self.propiedades['ac']
        cont_c = np.zeros(len(self.p_profundidad))
        for i, z in enumerate(self.p_profundidad):
            cont_c[i] = math.exp(-ac*z)
        self.cont_carbon = cont_c
        self.df_salida['fz'] = np.flip(cont_c)
        return cont_c
    
    
    def f_SwNit(self):
        WFP = self.WFP
        fs  = self.propiedades['fs']
        fwp = self.propiedades['fwp']
        swp = self.propiedades['swp']
        e2  = self.propiedades['e2']	 # e2
        e3  = self.propiedades['e3']
        sl  = self.propiedades['sl']
        sh  = self.propiedades['sh']
        
        fsw_nt = np.zeros(len(WFP))
        
        for i, wfp in enumerate(WFP):
            if sh < wfp <= 1:
                fsw_nt[i] = fs + (1-fs)*math.pow((1-wfp)/(1-sh),e2)
            elif sl <= wfp <= sh:
                fsw_nt[i] = 1
            elif swp <= wfp <= sl:
                fsw_nt[i] = fwp + (1-fwp)*math.pow((wfp-swp)/(sl-swp),e3)
        self.fsw_nt = fsw_nt
        self.df_salida['fsw_nt'] = np.flip(fsw_nt)
        return fsw_nt
        
      
    def f_SwDNit(self):
        WFP  = self.WFP
        sdn  = self.propiedades['sdn']
        ednt = self.propiedades['ednt']
        
        fsw_dnt = np.zeros(len(WFP))
        
        for i, wfp in enumerate(WFP):
            if wfp < sdn:
                fsw_dnt[i] = 0
            else:
                fsw_dnt[i] = math.pow((wfp-sdn)/(1-sdn),ednt)
        self.fsw_dnt = fsw_dnt
        self.df_salida['fsw_dnt'] = np.flip(fsw_dnt)
        return fsw_dnt
    
 ##### Parámetros que se utilizan en el calculo de C   
    def f_retardacion(self):
        kd = self.propiedades['kd']
        tetha = self.p_humedad
        rho =  self.propiedades['rho']
        r = np.zeros(len(tetha))
        
        for i, h in enumerate(tetha):
            r[i] = 1 + (rho*kd/(h))
        self.R = r
        self.df_salida['R'] = np.flip(r)
        return r
    
    def fNITmax(self):
        kr_max = self.propiedades['kr_max']
        kr = self.fsw_nt * kr_max
        self.krmax = kr
        self.df_salida['krmax'] = np.flip(kr)
        return kr
    
    
    def fDNTmax(self):
        Vmax = self.propiedades['Vmax']
        Vm = self.fsw_dnt * Vmax
        self.Vmax = Vm
        self.df_salida['Vmax'] = np.flip(Vm)
        return Vm
    
    def run(self):
            #cálculo de profundidad
        self.perfil_z()
        #### HIDRAULICA Y CONTENIDO DE HUMEDAD
        #Cálculo de presiones
        self.calc_presion()# = 
        #perfil de humedad
        self.calc_humedad()
        #WFP
        self.fWFP()
        #FACTORES PARA LAS REACCIONES
        ## efecto del contenido de carbono
        self.fcont_c()
        ## efecto de sw en la nit
        self.f_SwNit()
        ## efecto de sw en la desn
        self.f_SwDNit()
        #ELEMENTOS PARA CALCULAR LAS CONCENTRACIONES
        ### factor de retardación
        self.f_retardacion()
        ### kr max
        self.fNITmax()
        ##Vmax
        self.fDNTmax()
        ## Concentraciones  

    
    def get_data(self):
        return self.df_salida
    
    
    
    
class Estacionario(Modelo):
    def __init__(self, D=100, dz=1, propiedades=None, NH4=60, NO3=0.1):
        super().__init__(D, dz, propiedades, NH4, NO3)
        self.Td = self.incremento/self.vz
        self.mass_fl = None
        self.C_NH4 = None
        self.C_NO3 = None
        self.Ctotal = None
        self.allow_graph = False
        ## Concentraciones 
        
    
    def fconc(self, km, Co, mu_max, R=1, fz=1):
        td = self.Td
        ## Newton rhapson
        ini = 0.00001
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
    

    def Cnit(self):
        #parametros
        C_ini_NH4 = self.NH4
        C_ini_NO3 = self.NO3
        R_nt = self.df_salida['R']
        krmax = self.df_salida['krmax']
        vmax = self.df_salida['Vmax']       
        
        n = np.shape(self.p_profundidad)[0] 
        c_NH4 = np.zeros(n)
        c_NH4[0] = C_ini_NH4
        c_NO3 = np.zeros(n)
        c_NO3[0] = C_ini_NO3
        
        for i in range(1,n):
            if C_ini_NH4 > 0.001:
                c_NH4[i] = self.fconc(self.propiedades['Km_nit'], C_ini_NH4, krmax[i], R_nt[i])
            else:
                c_NH4[i] = 0
            
            C_ini_NH4 = c_NH4[i]
            
            C_ini_NO3 = c_NO3[i-1] + c_NH4[i-1] - c_NH4[i]
            if C_ini_NO3 > 0:
    
                c_NO3[i] = self.fconc(self.propiedades['Km_dnt'],C_ini_NO3, vmax[i])
            else:
                c_NO3[i] = 0.001
            
        self.C_NH4 = c_NH4
        self.df_salida['cNH4'] = c_NH4
        self.C_NO3 = c_NO3
        self.df_salida['cNO3'] = c_NO3
        self.Ctotal = c_NH4 + c_NO3    
        self.df_salida['Total N'] = self.Ctotal
        
        #Flujo de masa
        Ctotal = self.Ctotal
        hlr = self.propiedades['HLR']
        mass_fl = np.zeros(Ctotal.shape[0])
        for i in range(Ctotal.shape[0]):
            mass_fl[i] = (2/3) * Ctotal[i]* (2/3) * hlr / 100
        
        self.__mass_fl = mass_fl
            
            
    def ejecutar(self):
        self.run()
        ## Concentraciones
        self.Cnit()
        
        self.allow_graph = True
            
    def graficar(self, mass = False):
        if self.allow_graph:
            
            x = self.p_profundidad
            y1 = self.C_NH4
            y2 = self.C_NO3
            y3 = self.Ctotal
           
            
            
            fig, ax = plt.subplots()
            ax.plot(x, y1)
            ax.plot(x, y2)
            ax.plot(x, y3)
            ax.set_xlabel('Profundidad')
            ax.set_ylabel('Concentración')
            ax.set_title('Concentración N')
            ax.legend(['CNH4','CNO3','CTotal'])
            
            if mass:
                y_fl = self.mass_fl
                fig1, ax1 = plt.subplots()
                ax1.plot(x, y_fl)
                ax1.set_xlabel('Profundidad')
                ax1.set_ylabel('Mass flux $gm^2 d^{-1}$')
                ax1.set_title('Mass Flux TN')
                ax1.legend(['Mass Flux TN'])
                plt.show()
        else:
            print('Se debe ejecutar el modelo')
            
class Transitorio(Modelo):
    
    def __init__(self, Depth=100, dz=1, step = 10, dt = 1, propiedades=None, 
                 NH4=60, NO3=0.1, dispersion=0.001, reac=False, retardacion=True):
        super().__init__(Depth, dz, propiedades, NH4, NO3)   
        self.dt = dt
        self.pasos = step
        self.allow_graph = False
        self.dispersion = dispersion
        self.soluciones = None
        self.soluciones_NO3 = None 
        self.allow_reac = reac
        self.retardacion = retardacion
    
        
    def ejecutar(self):
        self.run()
        ##Definición de parametros
        L = self.p_profundidad.shape[0]
        p_retardacion = self.df_salida['R']
        reac_nit = self.df_salida['krmax']
        v_flujo = self.propiedades['HLR']
        km = self.propiedades['Km_nit'] 
        ##Parámetros
        dt = self.dt
        dz = self.incremento
        T = self.pasos
        D = self.dispersion
        ### Esquema temporal
        eps = 0.5
        
        
        #upwind
        peclet = (v_flujo*dz)/(2*D)
        print(peclet)
        lamb = 1 -0.78#(1/peclet) #+ (2/(np.e**(2*peclet)-1)) 
        ###Difusión artificial
        Dh = D*(1 + lamb*peclet)
        #condiciones iniciales y de frontera
        condicion = np.zeros(L)
        frontera = self.NH4
        condicion[0] = frontera
        
        incognitas = L-1
        
        react_array = np.zeros(incognitas)
        ###
        
        soluciones = [np.array([condicion])] 

        for t in range(T):#ciclo de tiempo
        
            matriz = np.zeros((incognitas, incognitas))
            vector = np.zeros(incognitas)
            
            for x in range(incognitas):
                if self.allow_reac:
                    f_reac = lambda con, mu_max, km: (mu_max*con)/(km+con)
                    #Reacciones nitrificacion
                    reac_nt = f_reac(condicion[x+1], reac_nit[x+1], km)
                    if reac_nt > 0.0001:
                        react_array[x] = reac_nt
                    else:
                        react_array[x] = 0 
                    
                if self.retardacion:
                    
                    R = p_retardacion[x+1]
                else:
                    R = 1
                
                
                alfa = (Dh * dt)/(R*dz**2)
                beta = (v_flujo * dt)/(R*2 * dz)
                
                
                
                ## Llenado de la matriz, y el vector de coeficientes
                
                if x == 0:
                    matriz[x, x] = 1 + 2*eps*alfa
                    matriz[x, x+1] = -eps*(alfa - beta)
                    
                    nodo_b = condicion[x]*(1-eps)*(alfa+beta)
                    nodo_c = condicion[x+1]*(1-2*(1-eps)*alfa)
                    nodo_f = condicion[x+2]*(1-eps)*(alfa-beta)
                    vector[x] =  nodo_b + nodo_c + nodo_f + eps*(alfa + beta)*frontera# -reac_nt

                
                elif 0 < x < incognitas-1:
                    matriz[x, x+1] = -eps*(alfa - beta) 
                    matriz[x, x] = 1 + 2*eps*alfa
                    matriz[x, x-1] = -eps*(alfa + beta) 
                    
                    nodo_b = condicion[x]*(1-eps)*(alfa+beta)
                    nodo_c = condicion[x+1]*(1-2*(1-eps)*alfa)
                    nodo_f = condicion[x+2]*(1-eps)*(alfa-beta)                    
                    vector[x] =  nodo_b + nodo_c + nodo_f #-reac_nt
                    
                    

                elif x == incognitas-1:
                    matriz[x, x-1] = -eps*(alfa + beta)
                    matriz[x, x] = 1 + eps*(alfa + beta)
                    
                    nodo_b = condicion[x]*(1-eps)*(alfa+beta)
                    nodo_c = condicion[x+1]*(1-(1-eps)*(alfa+beta))                    
                    
                    
                    vector[x] = nodo_b + nodo_c #-reac_nt
                    

                    
            ##Agregar algoritmo de gradiente biconjugado
            
            sol = np.linalg.solve(matriz, vector-react_array) #- react_array #

            
            
            condicion[1:] = sol[:]
            soluciones.append(np.array([condicion]))
        
        self.soluciones = soluciones 
        
        
        
        
    def graficar(self):
        
        plt.plot(self.p_profundidad, self.soluciones[-1].T, label='Upwind')
        #plt.plot(self.p_profundidad, self.soluciones_NO3.T)


if __name__ == '__main__':
    #Diccionario de propiedades
    dic_h = {
    'HLR':0.245,
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
    
    'kr_max':0.010,#	  kr max
    'Km_nit':0.50,#3	  Km,nit
    'bnit':0.35,#	  bnit 
    'sl':0.67,#	 sl
    'sh':0.81,#	  sh
    ####### dnt
    'ednt':3.77,#	  ednt
    'Vmax':0.056,#	  Vmax 
    'Km_dnt':50.00,#	  Km,dnt
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
    dat = mod.get_data()
    
    #est = Estacionario(120, 1, propiedades=dic_h, NH4=10, NO3 = 1)
    #est.ejecutar()
    #e_dat = est.get_data()
    #est.graficar()
    
    tran = Transitorio(1000,1,2000,1,dic_h,
                       #reac=True,
                       retardacion=False,
                       NH4=1,
                       dispersion=0.005
                       )
    tran.ejecutar()
    tran.graficar()
    
    
    
