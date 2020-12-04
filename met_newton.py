# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 14:26:01 2019

@author: Jorge Antonio Matías López
"""

#Newton
import math

def newton(ini):
    con = lambda c : c + math.log(c) + 10 
    der_con = lambda c: 1 + 1/c
    it = 1
    while True:
        nv = ini - con(ini)/der_con(ini)

        if abs(ini - nv) < 0.0000001 or it>50:
            return nv, it
            break
        else:
            ini = nv
            it += 1

res = newton(0.00001)