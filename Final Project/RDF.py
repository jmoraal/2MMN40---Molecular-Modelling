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
    INPUT: .xyz file, desired number of timesteps and indicator of how many steps to skip inbetween
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

def plotHist(*data, sizexAxis, nrBins):
    global binSize, bins, ballVolumes, binVolumes, histo, normalizedHisto
    fig = plt.figure(figsize =(10, 7)) 
    binSize = sizexAxis/nrBins
    bins = np.arange(0, sizexAxis, binSize)
    ballVolumes = (4/3)*np.pi*bins**3
    binVolumes = (np.roll(ballVolumes, -1) - ballVolumes)[0:len(bins)-1]
    
    # print(bins)
    # print(binVolumes)
    # print(binVolumes*rho)
    
    for i in range(0, len(data)):
        histo = np.histogram(data[i], bins = bins, density=False)
        # print(histo)
        normalizedHisto = histo[0]/(binVolumes*rho*nrOfMolecules*nrOfTimeSteps)  #nrOfMolecules) # normalize for volume of shell and average density. Also, avg over number of timesteps
        # print(normalizedHisto)
        labels = ["O-O", "O-H"]
        plt.plot(bins[0:len(bins)-1], normalizedHisto, label=labels[i])
    plt.xlabel("r $(\AA)$")
    plt.ylabel("g(r)")
    plt.legend(prop={'size': 16})
    
def RDFPerTimeStep(x, timeStep):
    global OwaterInd, HwaterInd, rOwaterHwater
    dist = distAtomsPBC(x[i,:,:])
    #dist[np.where(notInSameMolecule == False)] = 0
    
    OwaterInd = np.array(np.where((types == 0) & (molecule == 0))).flatten()
    HwaterInd = np.array(np.where((types == 1) & (molecule == 0))).flatten()
    OethanolInd = np.array(np.where((types == 0) & (molecule == 1))).flatten()
    HethanolInd = np.array(np.where((types == 1) & (molecule == 1))).flatten()
    
    rOwater = dist[OwaterInd,:]
    rOethanol = dist[OethanolInd,:]
    
    # rOwaterOwater = np.triu(rOwater[:,OwaterInd]) # same indices means double count?? 
    # rOethanolOethanol = np.triu(rOethanol[:,OethanolInd])
    rOwaterOwater = rOwater[:,OwaterInd] 
    rOethanolOethanol = rOethanol[:,OethanolInd]
    rOethanolOwater = rOethanol[:,OwaterInd]
    rOwaterHwater = rOwater[:,HwaterInd]
    rOethanolHethanol = rOethanol[:,HethanolInd]
    rOethanolHwater = rOethanol[:,HwaterInd]
    
    rOwaterOwater = rOwaterOwater[np.nonzero(rOwaterOwater)]
    rOethanolOethanol = rOethanolOethanol[np.nonzero(rOethanolOethanol)]
    rOethanolOwater = rOethanolOwater[np.nonzero(rOethanolOwater)]
    rOwaterHwater = rOwaterHwater[np.nonzero(rOwaterHwater)]
    rOethanolHethanol = rOethanolHethanol[np.nonzero(rOethanolHethanol)]
    rOethanolHwater = rOethanolHwater[np.nonzero(rOethanolHwater)]
    
    return(rOwaterOwater, rOethanolOethanol, rOethanolOwater, rOwaterHwater, rOethanolHethanol, rOethanolHwater)


outputFileName = "Water31.08ThermostatOutput.xyz"
topologyFileName = "Water31.08Topology.txt"
distAtomsPBC.boxSize = 31.08
rho = 0.032592 # molecules per Anstrom^3 (0.0999 atoms per Angstrom^3)
nrOfTimeSteps = 20

# outputFileName = "Ethanol32.22ThermostatOutput.xyz"
# topologyFileName = "Ethanol32.22Topology.txt"
# distAtomsPBC.boxSize = 32.22
# rho = 0.0104896*10**4 # particles per Anstrom^3
# nrOfTimeSteps = 10

# outputFileName = "MixedMoleculesOutput.xyz"
# topologyFileName = "MixedMoleculesTopology.txt"
# distAtomsPBC.boxSize = 31.08
# rho = 0.032592 # particles per Anstrom^3
# nrOfTimeSteps = 1


# outputFileName = "Water150Output.xyz"
# topologyFileName = "Water150Topology.txt"
# distAtomsPBC.boxSize = 20
# rho = 0.032592 # particles per Anstrom^3 3.345
# nrOfTimeSteps = 10

types, x = readXYZOutput(outputFileName, nrOfTimeSteps, Nskip = 4)
nrOfMolecules, molecule, notInSameMolecule = readTopologyFile(topologyFileName)

OwOw = []
OeOe = [] 
OeOw = [] 
OwHw = [] 
OeHe = [] 
OeHw = []
for i in range(0,nrOfTimeSteps):
    rOwaterOwater, rOethanolOethanol, rOethanolOwater, rOwaterHwater, rOethanolHethanol, rOethanolHwater = RDFPerTimeStep(x, i)
    OwOw.append(rOwaterOwater.flatten())
    OeOe.append(rOethanolOethanol.flatten())
    OeOw.append(rOethanolOwater.flatten())
    OwHw.append(rOwaterHwater.flatten())
    OeHe.append(rOethanolHethanol.flatten())
    OeHw.append(rOethanolHwater.flatten())

nrBins = 100
sizexAxis = 10

plotHist(OwOw, OwHw, sizexAxis = sizexAxis, nrBins = nrBins)

# plotHist(OeOe, OeHe, sizexAxis = sizexAxis, nrBins = nrBins)

    
    
    
### Hieronder de vorige versie, aub nog even niet wissen ###
    
# def distAtomsPBC(x, x2):
#     """ Computes distances between all atoms in closest copies, taking boundaries into account"""   
#     global dist, diff
#     diff = x - x2[:,np.newaxis]
#     diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 
#     dist = np.linalg.norm(diff,axis = 2)
    
#     return(dist.tolist()) 
# def plotDensity(*data):
#     fig = plt.figure(figsize =(10, 7)) 
#     labels = ["O-O", "O-H"]
#     for i in range(0, len(data)):
#         sns.distplot(data[i], hist = False, kde = True, label = labels[i])
#     plt.legend(prop={'size': 16})
#     plt.title("Radial Distribution Density")
#     plt.xlabel("r $(\AA)$")
#     plt.ylabel("Density")
    
    
# outputFileName = "Water31.08ThermostatOutput.xyz"
# topologyFileName = "Water31.08Topology.txt"
# distAtomsPBC.boxSize = 31.08
# rho = 0.32592 # particles per Anstrom^3

# nrOfTimeSteps = 3
# # timeStep = 1 # in range [0,nrOfTimeSteps - 1]

# types, x = readXYZOutput(outputFileName, nrOfTimeSteps)
# molecule, notInSameMolecule = readTopologyFile(topologyFileName)

# xOwater = x[:,np.where((types == 0) & (molecule == 0)),:]
# xOethanol = x[:,np.where((types == 0) & (molecule == 1)),:]
# xHwater = x[:,np.where((types == 1) & (molecule == 0)),:]
# xHethanol = x[:,np.where((types == 1) & (molecule == 1)),:]

# nrOwater = np.shape(xOwater)[2]
# nrOethanol = np.shape(xOethanol)[2]
# nrHwater = np.shape(xHwater)[2]
# nrHethanol = np.shape(xHethanol)[2]

# rOwaterOwater = []
# rOethanolOethanol = []
# rOwaterOethanol = []
# rOwaterHwater = []
# rOethanolHethanol = []
# rOethanolHwater = []

# for i in range(0,nrOfTimeSteps):
#     xOwaterAtTimei = xOwater[i,:,:,:].reshape(nrOwater,3)
#     xOethanolAtTimei = xOethanol[i,:,:,:].reshape(nrOethanol,3)
#     xHwaterAtTimei = xHwater[i,:,:,:].reshape(nrHwater,3)
#     xHethanolAtTimei = xHethanol[i,:,:,:].reshape(nrHethanol,3)
    
#     rOwaterOwater.append(distAtomsPBC(xOwaterAtTimei, xOwaterAtTimei))
#     rOethanolOethanol.append(distAtomsPBC(xOethanolAtTimei, xOethanolAtTimei))
#     rOwaterOethanol.append(distAtomsPBC(xOwaterAtTimei, xOethanolAtTimei))
#     rOwaterHwater.append(distAtomsPBC(xOwaterAtTimei, xHwaterAtTimei))
#     rOethanolHethanol.append(distAtomsPBC(xOethanolAtTimei, xHethanolAtTimei))
#     rOethanolHwater.append(distAtomsPBC(xOethanolAtTimei, xHwaterAtTimei))

# rOwaterOwater = np.asarray(rOwaterOwater).flatten()
# rOethanolOethanol = np.asarray(rOethanolOethanol).flatten()
# rOwaterOethanol = np.asarray(rOwaterOethanol).flatten()
# rOwaterHwater = np.asarray(rOwaterHwater).flatten()
# rOethanolHethanol = np.asarray(rOethanolHethanol).flatten()
# rOethanolHwater = np.asarray(rOethanolHwater).flatten()


# plotHist(rOwaterOwater, binSize)
# plotHist(rHwaterHawter, binSize)
# nrBins = 200
# sizexAxis = 6
# labels = ["O-O", "O-H"]
# plotHist(rOwaterOwater, rOwaterHwater, sizexAxis = sizexAxis, nrBins = nrBins, labels = labels)
# plotHist(rOwaterOwater, sizexAxis = sizexAxis, nrBins = nrBins)

# plt.hist(rOwaterHwater, bins = 'auto')

# print(np.max(rOwaterOwater) <= np.sqrt(2*(distAtomsPBC.boxSize/2)**2 + (distAtomsPBC.boxSize/2)**2)) 

# plotDensity(rOwaterOwater, rOwaterHwater, rHwaterHwater)
# plotDensity(rHwaterHwater)





























