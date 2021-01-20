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
    
        molecule = np.zeros(len(types))
        atom = 0
        for i in range(1,nrOfMolecules+1):
            for j in range(0,len(lines[i].split())):
                if len(lines[i].split()) == 3:
                    molecule[atom] = 0
                elif len(lines[i].split()) == 9:
                    molecule[atom] = 1
                atom += 1       
        
        # molecule = molecule.tolist()
        # molecule = list(map(int, molecule))
        return(molecule)
    
def readXYZOutput(fileName): 
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
    
    nrOfTimeSteps = 1
    
    atomTypesList = firstColumn[2:2+nrOfAtoms]
    typesDict = {'O': 0, 'H': 1, 'C': 2} # 0 oxygen, 1 hydrogen, 2 carbon
    atomTypes = np.vectorize(typesDict.get)(atomTypesList)
    
    atomPositions = np.zeros([nrOfTimeSteps,nrOfAtoms,3])
    
    lineNr = len(lines) - nrOfTimeSteps*(nrOfAtoms + 2)# start reading here, until the end
    
    i = 0
    while lineNr < len(lines):
        atomPositions[i,:,:] = np.asarray(lines[lineNr+2:lineNr+2+nrOfAtoms], dtype=float)
        lineNr = lineNr + nrOfAtoms + 2
        i += 1
        
    return(atomTypes, atomPositions)

def distAtomsPBC(x, x2):
    """ Computes distances between all atoms in closest copies, taking boundaries into account"""   
    diff = x - x2[:,np.newaxis]
    diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 
    dist = np.linalg.norm(diff,axis = 2)
    return(dist) 

outputFileName = "Water29Output.xyz"
topologyFileName = "Water29Topology.txt"

# outputFileName = "MixedMoleculesOutput.xyz"
# topologyFileName = "MixedMoleculesTopology.txt"

outputFileName = "Mixture29Output.xyz"
topologyFileName = "Mixture29Topology.txt"

distAtomsPBC.boxSize = 29
types, x = readXYZOutput(outputFileName)
molecule = readTopologyFile(topologyFileName)

xOwater = x[:,np.where((types == 0) & (molecule == 0)),:]
xOethanol = x[:,np.where((types == 0) & (molecule == 1)),:]
xHwater = x[:,np.where((types == 1) & (molecule == 0)),:]
xHethanol = x[:,np.where((types == 1) & (molecule == 1)),:]

nrOwater = np.shape(xOwater)[2]
nrOethanol = np.shape(xOethanol)[2]
nrHwater = np.shape(xHwater)[2]
nrHethanol = np.shape(xHethanol)[2]

timeStep = 0
xOwater = xOwater[timeStep,:,:,:].reshape(nrOwater,3)
rOwaterOwater = distAtomsPBC(xOwater, xOwater)
xOethanol = xOethanol[timeStep,:,:,:].reshape(nrOethanol,3)
rOwaterOethanol = distAtomsPBC(xOwater, xOethanol)
xHwater = xHwater[timeStep,:,:,:].reshape(nrHwater,3)
rHwaterOethanol = distAtomsPBC(xHwater, xOethanol)

rOwaterOwater = np.around(rOwaterOwater,2).reshape(nrOwater*nrOwater)
rOwaterOethanol = np.around(rOwaterOethanol,2).reshape(nrOwater*nrOethanol)
rHwaterOethanol = np.around(rHwaterOethanol,2).reshape(nrHwater*nrOethanol)

dr = distAtomsPBC.boxSize/10

fig = plt.figure(figsize =(10, 7)) 
bins = [0, dr, 2*dr, 3*dr, 4*dr, 5*dr, 6*dr, 7*dr, 8*dr, 9*dr, 10*dr]
plt.hist(rOwaterOwater, bins = bins)

plt.title("RDF") 
plt.show()


fig = plt.figure(figsize =(10, 7)) 
plt.hist(rHwaterOethanol, bins = bins, color = "black")
# TODO: shapes of the two histograms are too identical, cant be right


















