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
    return(kinetic,  bondPot, anglePot, dihedralPot, LJpot, tempBefore, tempAfter)
    
kinetice, bondPote, anglePote, dihedralPote, LJpote, tempBeforee, tempAftere = readMeasurables("Ethanol32.22ThermostatMeasurables.txt")
kineticw, bondPotw, anglePotw, dihedralPotw, LJpotw, tempBeforew, tempAfterw = readMeasurables("Water31.08ThermostatMeasurables.txt")
kineticm, bondPotm, anglePotm, dihedralPotm, LJpotm, tempBeforem, tempAfterm = readMeasurables("Mixture32.29ThermostatMeasurables.txt")

Epotw = bondPotw+ anglePotw+ dihedralPotw+ LJpotw
Etotw = Epotw + kineticw
Epote = bondPote+ anglePote+ dihedralPote+ LJpote
Etote = Epote + kinetice
Epotm = bondPotm+ anglePotm+ dihedralPotm+ LJpotm
Etotm = Epotm + kineticm



def summary(data):
    data = data[-int(len(data)/10):] # only looks at last 10% of data to avoid warmup period
    # print('Avg: ',np.average(data))
    # print('St. dev.: ',np.std(data))
    # print('Max: ', np.max(data))
    # print('Min: ', np.min(data))
    return np.average(data)

# wat = summary(Epotw)/summary(Etotw)
# eth = summary(Epote)/summary(Etote)
# print(wat, eth)
# print(summary(Etote) / (32.22**3))
# print(summary(Etotw) / (31.08**3))
# print(summary(Etotm) / (32.29**3))
print(summary(Etote) / (3089))
print(summary(Etotw) / (3000))
print(summary(Etotm) / (939*3 + 61*9))


### PLOT ENERGY ###
def plotEnergy(fileName, dt): 
    measurables = readMeasurables(fileName)
    Ekin = measurables[0]
    Epot = np.transpose(np.array(measurables[1:4]))
    EpotTotal = np.sum(Epot, axis = 1)
    plt.clf() # Clears current figure
    plt.rcParams.update({'font.size': 12})
    t = np.arange(0,len(Ekin)*dt, dt)
    Etot = np.array(Ekin) + np.array(EpotTotal)
    plt.plot(t, Ekin, label = 'Kinetic energy')
    plt.plot(t, EpotTotal, label = 'Potential energy')
    plt.plot(t, Etot, label = 'Total energy')
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (AMU Å² / ps²)')
    plt.legend()
    plt.show()
    plt.savefig(fileName + "Energies.pdf", bbox_inches = 'tight')
    


def plotPotentials(fileName, dt): 
    global Epot
    measurables = readMeasurables(fileName)
    Epot = np.transpose(np.array(measurables[1:5]))
    #EpotTotal = np.sum(Epot, axis = 1)
    plt.clf() # Clears current figure
    plt.rcParams.update({'font.size': 12})
    t = np.arange(0,len(Epot)*dt, dt)
    # plt.plot(t, Epot)
    plt.plot(t, Epot[:,0], 'C0', t, Epot[:,1], 'C1', t, Epot[:,3], 'C3') #for water, which does not have dihedral potential
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (AMU Å² / ps²)')
    plt.legend( ('Bonds', 'Angles', 'Lennard-Jones'), loc = 'best') #for water
    # plt.legend( ('Bonds', 'Angles', 'Dihedrals', 'Lennard-Jones'), loc = 'best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    plt.show()
    plt.savefig(fileName + "Potentials.pdf", bbox_inches = 'tight')
    


# plotEnergy("Ethanol32.22" + "ThermostatMeasurables.txt", 0.003)
# plotPotentials("Ethanol32.22" + "ThermostatMeasurables.txt", 0.003)

# plotEnergy("Mixture32.29" + "ThermostatMeasurables.txt", 0.002)
# plotPotentials("Mixture32.29" + "ThermostatMeasurables.txt", 0.002)

# plotEnergy("Water31.08" + "ThermostatMeasurables.txt", 0.004)
# plotPotentials("Water31.08" + "ThermostatMeasurables.txt", 0.004)

# plotEnergy("Ethanol_ex_noThermMeasurables.txt", 0.002)
# plotEnergy("Ethanol_ex_thermMeasurables.txt", 0.002)