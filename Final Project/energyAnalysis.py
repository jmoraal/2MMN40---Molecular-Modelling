# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 21:11:42 2021

@author: s161981
"""

import numpy as np
from matplotlib import pyplot as plt 

def readMeasurables(fileName): 
    global lines
    """Read a measurables file. (Les Measurables, haha)"""
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
            
    #IMPORTANT: for older runs, return MINUS dihedralPot! Original code contained wrong sign
    return(kinetic,  bondPot, anglePot, -dihedralPot, LJpot, tempBefore, tempAfter)
    
#kinetic, bondPot, anglePot, dihedralPot, LJpot, tempBefore, tempAfter = readMeasurables("Ethanol32.22ThermostatMeasurables.txt")


### PLOT ENERGY ###
def plotEnergy(fileName, dt): 
    measurables = readMeasurables(fileName)
    Ekin = measurables[0]
    Epot = np.transpose(np.array(measurables[1:4]))
    EpotTotal = np.sum(Epot, axis = 1)
    plt.clf() # Clears current figure
    t = np.arange(0,len(Ekin)*dt, dt)
    Etot = np.array(Ekin) + np.array(EpotTotal)
    plt.plot(t, Ekin, label = 'Kinetic energy')
    plt.plot(t, EpotTotal, label = 'Potential energy')
    plt.plot(t, Etot, label = 'Total energy')
    plt.title('Energy in the system')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (AMU Å² / ps²)')
    plt.legend()
    plt.show()
    plt.savefig(fileName + "Energies.pdf")
    
plotEnergy("Ethanol32.22ThermostatMeasurables.txt", 0.003)


def plotPotentials(fileName, dt): 
    global Epot
    measurables = readMeasurables(fileName)
    Epot = np.transpose(np.array(measurables[1:5]))
    #EpotTotal = np.sum(Epot, axis = 1)
    plt.clf() # Clears current figure
    t = np.arange(0,len(Epot)*dt, dt)
    plt.plot(t, Epot)
    plt.title('Levels of potential energy')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (AMU Å² / ps²)')
    plt.legend( ('Bonds', 'Angles', 'Dihedrals', 'Lennard-Jones'), loc = 'best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    plt.show()
    plt.savefig(fileName + "Potentials.pdf")
    
plotPotentials("Ethanol32.22ThermostatMeasurables.txt", 0.003)
