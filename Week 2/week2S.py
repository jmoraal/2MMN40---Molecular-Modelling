# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981

NOTE: using tab for indentation
"""

import numpy as np

### WEEK 1 ###
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

#testDist = distAtTime(testPos)

### WEEK 2 ###

# TODO
# - in de instructie staat dat er een betere oplossing is voor dot product in FAngleOnAtoms
#    Update: alternatieve uitdrukking voor t, is die beter?
# - waarschijnlijk zijn er efficientere methodes voor FTotalOnAtoms en dat ding printen
# - wat Ruben zei in Teams over volgorde van atomen (in group X donderdag 19 nov)

# BOND
def Vbond(r):
    return(1/2*k*(r-r0)**2)

def Fbond(r):
    return(-k*(r-r0))

def FBondOnAtoms(a,b):
    r = np.linalg.norm(a-b)
    Fa = Fbond(r)*(a-b)/np.linalg.norm(a-b)
    return(np.asarray([Fa, -Fa]))

# ANGLE
def Vangle(t):
    return 1/2*kt*(t-t0)**2

def Fangle(t):
    return -kt*(t-t0)

def FAngleOnAtoms(a,b,c):
    """Compute angular forces on 3-body atom.
    
    Mind the order of the arguments. The middle argument is supposed to be the middle atom (O in water).
    INPUT: positions of atoms a,b,c
    OUTPUT: angular force acting on each of the atoms
    """
    t = np.arccos(np.dot((a-b),(c-b))/(np.linalg.norm(a-b)*np.linalg.norm(c-b)))
    """ Alternative computation of t using cosine rule:
    ab = np.linalg.norm(a-b)
    bc = np.linalg.norm(c-b)
    ac = np.linalg.norm(a-c)
    t = np.arccos((ab**2 + bc**2 - ac**2)/(2*ab*bc))
    """
    normalVecA = np.cross(a-b,np.cross(a-b,c-b))
    normalVecC = np.cross(b-c,np.cross(a-b,c-b))
    Fa = Fangle(t)/np.linalg.norm(a-b) * normalVecA/np.linalg.norm(normalVecA)
    Fc = Fangle(t)/np.linalg.norm(c-b) * normalVecC/np.linalg.norm(normalVecC)
    return(np.asarray([Fa, -Fa-Fc, Fc]))

def FTotalOnAtoms(a,b,c):
    """Compute total forces on 3-body atom."""
    FBondAB = FBondOnAtoms(a,b)
    FBondBC = FBondOnAtoms(b,c)
    FAngle = FAngleOnAtoms(a,b,c)
    Fa = FBondAB[0] + FAngle[0]
    Fc = FBondBC[1] + FAngle[2]
    Fb = FBondAB[1] + FBondBC[0] + FAngle[1]
    return(np.asarray([Fb, Fa, Fc]))

# hydrogen example
types, xyzs = readXYZfile("HydrogenSingle.xyz", 0)
k = 24531/(10**2) # in kJ / (mol A^2)
r0 = 0.74 # in Angstrom

print(FBondOnAtoms(xyzs[0],xyzs[1]))

# water example
types, xyzs = readXYZfile("WaterSingle.xyz", 0)
k = 502416/(10**2) # in kJ / (mol A^2)
r0 = 0.9572 # in Angstrom
kt = 628.02
t0 = np.deg2rad(104.52)

print(FTotalOnAtoms(xyzs[1] , xyzs[0], xyzs[2]))
