# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981

NOTE: using tab for indentation
"""

import numpy as np




def readXYZfile(fileName): 
    """Read a .xyz file.
    
    Reads entire file for arbitrary number of timesteps
    INPUT: .xyz file 
    OUTPUT:  atom types and array of xyz-coord for every atom at every timestep.
    """
    lines = []
    firstColumn = []

    with open(fileName, "r") as inputFile:
        for line in inputFile: 
            splittedLine = line.split()
            firstColumn.append(splittedLine[0])
            lines.append(splittedLine[1:4])
        
    nrOfAtoms = int(firstColumn[0])
    timesteps = int(len(lines)/(nrOfAtoms + 2))
    #This works because for every block of nrOfAtoms positions, there are two other lines
    
    atomPositions = []
    atomTypes = []
    
    for i in range(timesteps):
        atomTypes.append(firstColumn[(2+(2+nrOfAtoms)*i):((2+nrOfAtoms)*(i+1))])
        atomPositions.append(lines[(2+(2+nrOfAtoms)*i):((2+nrOfAtoms)*(i+1))])
        
    atomPositions = np.asarray(atomPositions).astype(np.float)
    return(atomTypes,atomPositions)


def distAtTime(positions,timestep):
    posAtTime = positions[timestep]
    diff = posAtTime - posAtTime[:,np.newaxis]
    dist = np.linalg.norm(diff,axis = 2)
    return(dist)


def distAtTimeFromFile(fileName,timestep):
    posAtTime = readXYZfile(fileName)[1][timestep]
    diff = posAtTime - posAtTime[:,np.newaxis]
    dist = np.linalg.norm(diff,axis = 2)
    return(dist)

#waterSmallPos = readXYZfile("waterSmall.xyz")[1]
#print(distAtTime(waterSmallPos,0))
print(distAtTimeFromFile("waterSmall.xyz",0))
'''
vragen:
    - Is atomtypes of nrOfAtoms ooit niet constant?
    - Comments worden genegeerd; moeten we <t = ...> kunnen lezen?
    - Liever functies van position of positionAtTime?
'''

'''
Test!
'''



