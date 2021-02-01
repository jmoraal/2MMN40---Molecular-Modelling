# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 21:11:42 2021

@author: s161981
"""

import numpy as np
from matplotlib import pyplot as plt 

def readMeasurables(fileName): 
    global lines
    """Read a measurables file. Or: Les Measurables?"""
    with open(fileName, "r") as inputFile:
        lines = inputFile.readlines()
        
        nrTimesteps = len(lines)-1
        
        kinetic = np.zeros(nrTimesteps)
        bondPot = np.zeros(nrTimesteps)
        anglePot = np.zeros(nrTimesteps)
        dihedralPot = np.zeros(nrTimesteps)
        LJpot = np.zeros(nrTimesteps)
        tempBefore = np.zeros(nrTimesteps)
        tempAfter = np.zeros(nrTimesteps)
        
        for i in range(0,nrTimesteps):
            (kinetic[i],  bondPot[i], anglePot[i], dihedralPot[i], LJpot[i], tempBefore[i], tempAfter[i]) = lines[i+1].split()
                 
    return(kinetic,  bondPot, anglePot, dihedralPot, LJpot, tempBefore, tempAfter)
    
kinetic, bondPot, anglePot, dihedralPot, LJpot, tempBefore, tempAfter = readMeasurables("Ethanol32.22ThermostatMeasurables.txt")