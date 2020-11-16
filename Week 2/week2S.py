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

#testDist = distAtTime(testPos)

# hydrogen example
types, xyzs = readXYZfile("hydrogenSmall.xyz", 1)

k = 24531/100
r0 = 0.74 
xyzs = xyzs

def Vbond(r):
    return(1/2*k*(r-r0)**2)

def Fbond(r):
    return(-k*(r-r0))

def FBondOnAtoms(a,b):
    r = np.linalg.norm(a-b)
    Fa = Fbond(r)*(a-b)/np.linalg.norm(a-b)
    return(np.asarray([Fa, -Fa]))

# print(FBondOnAtoms(xyzs[0],xyzs[1]))

# water example
types, xyzs = readXYZfile("WaterExampleWeek2.xyz", 0)

k = 502416/100
r0 = 0.9572
kt = 628.02
t0 = 104.52

def Vangle(t):
    return 1/2*kt*(t-t0)**2

def Fangle(t):
    return -kt*(t-t0)

def FAngleOnAtoms(a,b,c):
    t = np.dot((a-b),(b-c))/(np.linalg.norm(a-b)*np.linalg.norm(b-c))
    rAB = np.linalg.norm(a-b)
    rBC = np.linalg.norm(b-c)
    # print(a-b)
    # print(b-c)
    # print(t)
    # print(rAB)
    normalVecA = np.cross(b-a,np.cross(b-a,b-c))
    normalVecC = np.cross(c-b,np.cross(b-a,b-c))
    Fa = Fangle(t)*normalVecA/(np.linalg.norm(normalVecA)*rAB)
    Fc = Fangle(t)*normalVecC/(np.linalg.norm(normalVecC)*rBC)
    forces = np.asarray([Fa, -Fa-Fc, Fc])
    return(forces)

# print(FBondOnAtoms(xyzs[0],xyzs[1]))
# print("...")
print(FAngleOnAtoms(xyzs[0], xyzs[1], xyzs[2]))
# FAngleOnAtoms(xyzs[0], xyzs[1], xyzs[2])













