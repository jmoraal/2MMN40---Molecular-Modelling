# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981

NOTE: using tab for indentation
"""

import numpy as np

"""
To do: 
    Adapt readXYZfile: 
    - read 1st line of file separately to obtain array dimensions 
    - take timeStep as argument and only return one [nrOfAtoms, 3] array
    - reduce memory usage
"""    
def readXYZnrAtoms(fileName): 
    with open(fileName, "r") as inputFile:
        firstString = inputFile.readline().split()
        nrOfAtoms = int(firstString[0])
    return nrOfAtoms

def readXYZfile(fileName, timeStep): 
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
    #timesteps = int(len(lines)/(nrOfAtoms + 2))
    #This works because for every block of nrOfAtoms positions, there are two other lines
    
    atomTypes = firstColumn[(2+(2+nrOfAtoms)*timeStep):((2+nrOfAtoms)*(timeStep+1))]
    atomPositions = lines[(2+(2+nrOfAtoms)*timeStep):((2+nrOfAtoms)*(timeStep+1))]
        
    atomPositions = np.asarray(atomPositions).astype(np.float)
    return(atomTypes,atomPositions)


def distAtTime(positions):
    """ Computes distances between all atoms at given timestep """
    diff = positions - positions[:,np.newaxis]
    dist = np.linalg.norm(diff,axis = 2)
    return(dist)

#testPos = readXYZfile("Methane.xyz", 1)[1]
#testDist = distAtTime(testPos)

""" Also an option, but now still inefficient because it reads
    the entire file and then takes one slice of the resulting array
def distAtTimeFromFile(fileName,timestep):
    posAtTime = readXYZfile(fileName)[1][timestep]
    diff = posAtTime - posAtTime[:,np.newaxis]
    dist = np.linalg.norm(diff,axis = 2)
    return(dist)
"""




