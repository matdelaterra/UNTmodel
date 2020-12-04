# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 21:53:42 2020

@author: jor_a
"""
def n2_a_letras(n):
    n = str(int(n))
    cant = ''
    if len(n)==2:
        if int(n)>29:
            cant = cant + decenas[int(n[-2])]
            if int(n[-1]) !=0:
                cant += ' y '+ iniciales[int(n[-1])]      
        else:
            if int(n[-2])!=0:
                cant += especiales[int(n)]
            else:
                cant += iniciales[int(n[-1])] 
    else:
        cant += iniciales[int(n)]  
    return cant

def cientos2num(n):
    n = str(int(n))
    if len(n) <=3:
        cantidades = ''
        if len(n)>2:
            if n[0] == '1':
                if int(n[-2:]) == 0:
                    cantidades +='cien'
                else:
                    cantidades += 'ciento'
            else:
                cantidades += centenas[int(n[0])]
            cantidades += ' ' + n2_a_letras(n[-2:])
        else:
            cantidades += n2_a_letras(n)
    return cantidades

def miles2letra(n):
    n = str(int(n))
    numero = ''
    if len(n)>3:
        miles = cientos2num(int(n[:-3]))
        cientos = cientos2num(int(n[-3:]))
        numero += miles + ' mil ' + cientos
    else:
        numero += cientos2num(n)
    return numero

iniciales = {0: '',
1: 'un',
2: 'dos',
3: 'tres',
4: 'cuatro',
5: 'cinco',
6: 'seis',
7: 'siete',
8: 'ocho',
9: 'nueve'}

especiales = {
10:'diez',
11:'once',
12:'doce',
13:'trece',
14:'catorce',
15:'quince',
16:'dieciseis',
17:'diecisiete',
18:'dieciocho',
19:'diecinueve',
20:'veinte',
21:'veintiun',
22:'veintidos',
23:'veintitres',
24:'veinticuatro',
25:'veinticinco',
26:'veintiséis',
27:'veintisiete',
28:'veintiocho',
29:'veintinueve'}
decenas={
3:'treinta',
4:'cuarenta',
5:'cincuenta',
6:'sesenta',
7:'setenta',
8:'ochenta',
9:'noventa'}
centenas = {
1:'ciento',
2:'doscientos',
3:'trescientos',
4:'cuatrocientos',
5:'quinientos',
6:'seiscientos',
7:'setecientos',
8:'ochocientos',
9:'novecientos',
}