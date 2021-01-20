#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:04:55 2021

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)
"""
import numpy as np
from matplotlib import pyplot as plt 

def readTopologyFile(fileNameTopology): 
    """Read a topology file."""
    with open(fileNameTopology, "r") as inputFile:
        lines = inputFile.readlines()
        
        nrOfMolecules = int(lines[0].split()[1])
    
        molecules = []
        for i in range(1,nrOfMolecules+1):
            molecules.append(list(map(int,lines[i].split())))
        
        notInSameMolecule = np.ones((len(types), len(types)), dtype=bool)
        for i,mol in enumerate(molecules):
            for j,at in enumerate(mol):
                for k,at2 in enumerate(mol):
                    # print(at, at2)
                    notInSameMolecule[at,at2] = False
        
        return(molecules, notInSameMolecule)
    
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
    
    nrOfTimeSteps = 4
    
    atomTypesList = firstColumn[2:2+nrOfAtoms]
    typesDict = {'O': 0, 'H': 1, 'C': 2}
    atomTypes = np.vectorize(typesDict.get)(atomTypesList)
    atomPositions = np.zeros([nrOfTimeSteps,nrOfAtoms,3])
    
    lineNr = int((nrOfAtoms+2)*(np.ceil((endTime - nrOfTimeSteps*dt - startTime)/dt))) # start reading here, until the end
    i = 0
    while i < nrOfTimeSteps:
        atomPositions[i,:,:] = np.asarray(lines[lineNr+2:lineNr+2+nrOfAtoms], dtype=float)
        lineNr = lineNr + nrOfAtoms + 2
        i += 1
        
    return(atomTypes, atomPositions)

def distAtomsPBC(x):
    """ Computes distances between all atoms in closest copies, taking boundaries into account"""    
    diff = x - x[:,np.newaxis] 
    diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 

    # idea: if dist > 0.5*boxsize in some direction (x, y or z), then there is a closer copy. 
    # subtracting 0.5*boxsize in every direction where it is too large yields direction vector to closest neighbour
    # Still check correctness!
    dist = np.linalg.norm(diff,axis = 2)
    return(diff,dist) 

outputFileName = "Water29Output.xyz"
topologyFileName = "Water29Topology.txt"

startTime = 0
endTime = 0.012
dt = 0.003
distAtomsPBC.boxSize = 29
types, x = readXYZOutput(outputFileName, startTime, endTime, dt)
molecules, notInSameMolecule = readTopologyFile(topologyFileName)

rOwaterOwater = x[:,np.where(types == 0),:] # TODO only water molecules not ethanol
nrOwaterOwater = np.shape(rOwaterOwater)[2]

time0 = rOwaterOwater[0,:,:,:].reshape(nrOwaterOwater,3)

distance = np.linalg.norm(time0 - time0[:,np.newaxis], axis = 2)
distance = np.around(distance,2).reshape(nrOwaterOwater*nrOwaterOwater)

dr = distAtomsPBC.boxSize/10

fig = plt.figure(figsize =(10, 7)) 
bins = [0, dr, 2*dr, 3*dr, 4*dr, 5*dr, 6*dr, 7*dr, 8*dr, 9*dr, 10*dr]
plt.hist(distance, bins = bins)

plt.title("RDF") 
plt.show()


















