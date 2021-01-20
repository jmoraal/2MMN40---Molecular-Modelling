# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:34:12 2021

@author: s161981
"""


import numpy as np

### Box size tester ###
def findSize(approxSize, substance):
    sizes = np.arange(approxSize - 3.0, approxSize + 3.0, 0.01)
    global d
    d = []
    for boxSize in sizes: 
        if (substance == 'water'):
            nrOfMolecules = int(33.328*((boxSize/10)**3)) # water
        elif(substance == 'ethanol'): 
            nrOfMolecules = int(10.272*((boxSize/10)**3)) # ethanol
        elif(substance == 'mixture'):
            nrOfMolecules = int(29.728*((boxSize/10)**3)) # mixture
        else: 
            print('Substance not known')            
        
        n = int(nrOfMolecules**(1/3)) + 1
        if (n**3 - nrOfMolecules ==0): 
            d.append(boxSize)
    return d

print(findSize(30,'ethanol'))