#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:04:55 2021

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)
"""
import time as timer
import numpy as np

def readXYZOutput(fileName, startTime, endTime, dt): 
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
    # print(len(lines))
    # print(nrOfAtoms)
    
    nrOfTimeSteps = 3
    
    atomTypes = firstColumn[2:2+nrOfAtoms]
    atomPositions = np.zeros([nrOfTimeSteps,nrOfAtoms,3])
    
    lineNr = int((nrOfAtoms+2)*(np.ceil((endTime - nrOfTimeSteps*dt - startTime)/dt))) # start reading here, until the end
    
    i=0
    while i < nrOfTimeSteps:
        atomPositions[i,:,:] = np.asarray(lines[lineNr+2:lineNr+2+nrOfAtoms], dtype=float)
        lineNr = lineNr + nrOfAtoms + 2
        i += 1
        
    # atomPositions = np.asarray(atomPositions)#.astype(np.float) # .reshape(nrOfAtoms, nrOfTimeSteps)
    # print(atomPositions)
    return(atomTypes, atomPositions)



outputFileName = "MixedMoleculesOutput.xyz"
startTime = 0
endTime = 0.3 
dt = 0.003
types, x = readXYZOutput(outputFileName, startTime, endTime, dt)

