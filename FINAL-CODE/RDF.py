#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:04:55 2021

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)
"""
import numpy as np
from matplotlib import pyplot as plt 
import seaborn as sns

def readTopologyFile(fileNameTopology): 
    """Read a topology file."""
    with open(fileNameTopology, "r") as inputFile:
        lines = inputFile.readlines()
        
        nrOfMolecules = int(lines[0].split()[1])
    
        moleculeofAtom = np.zeros(len(types)) # 0 water, 1 ethanol
        atom = 0
        for i in range(1,nrOfMolecules+1):
            for j in range(0,len(lines[i].split())):
                if len(lines[i].split()) == 3:
                    moleculeofAtom[atom] = 0
                elif len(lines[i].split()) == 9:
                    moleculeofAtom[atom] = 1
                atom += 1       
        
        molecules = []
        for i in range(1,nrOfMolecules+1):
            molecules.append(list(map(int,lines[i].split())))
        
        # create boolean matrix to indicate which atoms are part of same molecule:
        # TODO: Probably can still be done more efficiently?
        notInSameMolecule = np.ones((len(types), len(types)), dtype=bool)
        for i,mol in enumerate(molecules):
            for j,at in enumerate(mol):
                for k,at2 in enumerate(mol):
                    notInSameMolecule[at,at2] = False
        
        return(nrOfMolecules, moleculeofAtom, notInSameMolecule)
    
def readXYZOutput(fileName, nrOfTimeSteps, Nskip = 1): 
    """Read a .xyz file.
    
    Reads rest of file starting at nrOfTimeSteps from the end
    INPUT: .xyz file, desired number of timesteps and (optionally) indicator of how many steps to skip inbetween
    OUTPUT:  atom types and array of xyz-coord for atoms for given number of timesteps
    """
    lines = []
    firstColumn = []

    with open(fileName, "r") as inputFile:
        for line in inputFile: 
            splittedLine = line.split()
            firstColumn.append(splittedLine[0])
            lines.append(splittedLine[1:4])
        
    nrOfAtoms = int(firstColumn[0])
    
    atomTypesList = firstColumn[2:2+nrOfAtoms]
    typesDict = {'O': 0, 'H': 1, 'C': 2} # 0 oxygen, 1 hydrogen, 2 carbon
    atomTypes = np.vectorize(typesDict.get)(atomTypesList)
    
    atomPositions = np.zeros([nrOfTimeSteps,nrOfAtoms,3])
    
    lineNr = len(lines) - Nskip * nrOfTimeSteps*(nrOfAtoms + 2)# start reading here, until the end
    
    #i = 0
    #while lineNr < len(lines):
    for i in range(0,nrOfTimeSteps): 
        atomPositions[i,:,:] = np.asarray(lines[lineNr+2:lineNr+2+nrOfAtoms], dtype=float)
        lineNr = lineNr + Nskip *(nrOfAtoms + 2)
        #i += 1
        
    return(atomTypes, atomPositions)


def distAtomsPBC(x):
    """ Computes distances between all atoms in closest copies, taking boundaries into account"""   
    global dist, diff
    diff = x - x[:,np.newaxis]
    diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 
    dist = np.linalg.norm(diff,axis = 2)
    
    return(dist) 

def plotHist(x, nrOfTimeSteps, sizexAxis, nrBins):
    '''Normalize and plot RDF histogram for given axis length and nr of bins 
    
    INPUT: array of all particle coordinates at different times, nr of timesteps to read, histogram axis size and nr of bins
    OUTPUT: None; plot saved to pdf
    '''
    plt.clf()
    OwOw = []
    OeOe = [] 
    OeOw = [] 
    OwHw = [] 
    OeHe = [] 
    OeHw = []
    for i in range(0,nrOfTimeSteps):
        print(i, " out of ", nrOfTimeSteps)
        rOwaterOwater, rOethanolOethanol, rOethanolOwater, rOwaterHwater, rOethanolHethanol, rOethanolHwater = RDFPerTimeStep(x, i, sizexAxis, nrBins)
        OwOw.append(rOwaterOwater.flatten())
        OeOe.append(rOethanolOethanol.flatten())
        OeOw.append(rOethanolOwater.flatten())
        OwHw.append(rOwaterHwater.flatten())
        OeHe.append(rOethanolHethanol.flatten())
        OeHw.append(rOethanolHwater.flatten())
    
    #fig = plt.figure(figsize =(10, 7)) 
    binSize = sizexAxis/nrBins
    bins = np.arange(0, sizexAxis, binSize)
    ballVolumes = (4/3.0)*np.pi*bins**3
    binVolumes = (np.roll(ballVolumes, -1) - ballVolumes)[0:len(bins)-1]
    
    # Alternatively, approximate with ball surface:
    # dr = bins[1] - bins[0]
    # binVolumes = 4*np.pi*bins[1:]**2 * dr
    
    # data = [OwOw,OwHw,OeOe,OeHe,OeOw,OeHw]
    data = [OeOe,OeHe]
    labels = ["O-O", "O-H"] #when not plotting mixture
    # labels = ["O-O water", "O-H water", "O-O ethanol", "O-H ethanol", "O-O ethanol-water", "O-H ethanol-water",]
    for i in range(len(data)):
        histo = np.sum(data[i], axis = 0)
        normalizedHisto = histo/(nrOfTimeSteps)
            
        # Volume/density compensation:
        nAvg = binVolumes * (len(types)/distAtomsPBC.boxSize**3) # Now compensates atom-specific density by additional factor in RDFPerTimeStep
        normalizedHisto = normalizedHisto/(nAvg)
        
        plt.plot(bins[0:len(bins)-1], normalizedHisto, label=labels[i])
    plt.xlabel("r $(\AA)$")
    plt.ylabel("g(r)")
    plt.legend()
    plt.savefig(name + "RDF.pdf", bbox_inches = 'tight')
    
    
def RDFPerTimeStep(x, timeStep, sizeAxis, nrBins):
    '''Analyzes set of coordinates at given timestep to prepare for histogram'''
    dist = distAtomsPBC(x[timeStep,:,:])
    dist[np.where(notInSameMolecule == False)] = 0
    
    OwaterInd = np.array(np.where((types == 0) & (molecule == 0))).flatten()
    HwaterInd = np.array(np.where((types == 1) & (molecule == 0))).flatten()
    OethanolInd = np.array(np.where((types == 0) & (molecule == 1))).flatten()
    HethanolInd = np.array(np.where((types == 1) & (molecule == 1))).flatten()
    
    rOwater = dist[OwaterInd,:]
    rOethanol = dist[OethanolInd,:]
    
    rOwaterOwater = rOwater[:,OwaterInd] 
    rOethanolOethanol = rOethanol[:,OethanolInd] 
    rOethanolOwater = rOethanol[:,OwaterInd]
    rOwaterHwater = rOwater[:,HwaterInd] 
    rOethanolHethanol = rOethanol[:,HethanolInd]
    rOethanolHwater = rOethanol[:,HwaterInd]
        
    binSize = sizexAxis/nrBins
    bins = np.arange(0, sizexAxis, binSize)
    
    #Normalize: compensate for nr of atoms at either end by dividing by atom-specific density. (Would be more at place in plotHist, but was easier implemented here)
    rOwaterOwaterCount = np.histogram(rOwaterOwater[np.nonzero(rOwaterOwater)], bins = bins, density=False)[0] * len(types) / len(OwaterInd)**2
    rOethanolOethanolCount = np.histogram(rOethanolOethanol[np.nonzero(rOethanolOethanol)], bins = bins, density=False)[0] * len(types) / len(OethanolInd)**2
    rOethanolOwaterCount = np.histogram(rOethanolOwater[np.nonzero(rOethanolOwater)], bins = bins, density=False)[0] * len(types) / (len(OethanolInd) * len(OwaterInd))
    rOwaterHwaterCount = np.histogram(rOwaterHwater[np.nonzero(rOwaterHwater)], bins = bins, density=False)[0] * len(types) / (len(OwaterInd)*len(HwaterInd)) 
    rOethanolHethanolCount = np.histogram(rOethanolHethanol[np.nonzero(rOethanolHethanol)], bins = bins, density=False)[0] * len(types) / (len(OethanolInd)*len(HethanolInd))  
    rOethanolHwaterCount = np.histogram(rOethanolHwater[np.nonzero(rOethanolHwater)], bins = bins, density=False)[0] * len(types) / (len(OethanolInd)*len(HwaterInd))   

    return(rOwaterOwaterCount, rOethanolOethanolCount, rOethanolOwaterCount, rOwaterHwaterCount, rOethanolHethanolCount, rOethanolHwaterCount)

### Cases to analyse - uncomment desired one ###
# name = "Water31.08Thermostat"
# topologyFileName = "Water31.08Topology.txt"
# distAtomsPBC.boxSize = 31.08
# nrOfTimeSteps = 100

# name = "Ethanol32.22Thermostat"
# topologyFileName = "Ethanol32.22Topology.txt"
# distAtomsPBC.boxSize = 32.22
# nrOfTimeSteps = 1000

name = "Mixture32.29Thermostat"
topologyFileName = "Mixture32.29Topology.txt"
distAtomsPBC.boxSize = 32.29
nrOfTimeSteps = 500

# name = "Water150"
# topologyFileName = "Water150Topology.txt"
# distAtomsPBC.boxSize = 20
# nrOfTimeSteps = 10


outputFileName = name + "Output.xyz"
types, x = readXYZOutput(outputFileName, nrOfTimeSteps)
nrOfMolecules, molecule, notInSameMolecule = readTopologyFile(topologyFileName)

nrBins = 100
sizexAxis = 10


plotHist(x,nrOfTimeSteps, sizexAxis, nrBins)
