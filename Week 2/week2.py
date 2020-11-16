# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 11:46:53 2020

@author: s161981
"""
import numpy as np
import week1 as xyz

def Vbond (r,r0,k):
    return 0.5*k*(r-r0)**2

def Fbond (r0,k):
    return -k*(-r0)

Fbond = Fbond(0.5,0.74,24531.0)

forceA = Fbond 