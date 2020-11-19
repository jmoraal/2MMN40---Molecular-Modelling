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

# BOND
def Vbond(r):
    return(1/2*k*(r-r0)**2)

def Fbond(r):
    return(-k*(r-r0))

def FBondOnAtoms(a,b):
    r = np.linalg.norm(b-a)
    Fa = Fbond(r)*(b-a)/np.linalg.norm(b-a)
    return(np.asarray([Fa, -Fa]))

# ANGLE
def Vangle(t):
    return 1/2*kt*(t-t0)**2

def Fangle(t):
    return -kt*(t-t0)

def FAngleOnAtoms(a,b,c):
    t = np.arccos(np.dot((a-b),(c-b))/(np.linalg.norm(a-b)*np.linalg.norm(c-b)))
    # print(np.rad2deg(t))
    rAB = np.linalg.norm(a-b)
    rBC = np.linalg.norm(c-b)
    # print(a-b)
    # print(b-c)
    # print(rAB)
    # print(rBC)
    # print(Fangle(t))
    normalVecA = np.cross(a-b,np.cross(a-b,c-b))
    # print(normalVecA)
    normalVecC = np.cross(c-b,np.cross(a-b,c-b))
    # print(np.cross(b-a,c-b))
    # print(np.dot(np.cross(b-a,c-b),c-b))
    # print(normalVecA)
    # print(normalVecC)
    Fa = Fangle(t)/rAB * normalVecA/np.linalg.norm(normalVecA)
    Fc = Fangle(t)/rBC * normalVecC/np.linalg.norm(normalVecC)
    forces = np.asarray([Fa, -Fa-Fc, Fc])
    return(forces)


# hydrogen example
types, xyzs = readXYZfile("HydrogenSingle.xyz", 0)
k = 24531/(10**2) # in kJ / (mol A^2)
r0 = 0.74 # in Angstrom

# print(FBondOnAtoms(xyzs[0],xyzs[1]))


# water example
types, xyzs = readXYZfile("WaterSingle.xyz", 0)
k = 502416/(10**2) # in kJ / (mol A^2)
r0 = 0.9572 # in Angstrom
kt = 628.02
t0 = np.deg2rad(104.52)

# print(FBondOnAtoms(xyzs[0],xyzs[1]))
print(FBondOnAtoms(xyzs[1], xyzs[0]))
print("...")
print(FAngleOnAtoms(xyzs[1], xyzs[0], xyzs[2])[0])
print("...")
print(FAngleOnAtoms(xyzs[1], xyzs[0], xyzs[2])[0] + FBondOnAtoms(xyzs[1], xyzs[0])[0])
# f = FAngleOnAtoms(xyzs[1], xyzs[0], xyzs[2]) # mind order

# print(FBondOnAtoms(xyzs[1], xyzs[0]))
# print(FBondOnAtoms(xyzs[0], xyzs[2]))

# print(FBondOnAtoms(xyzs[1], xyzs[0])[0] + FAngleOnAtoms(xyzs[1], xyzs[0], xyzs[2])[0])






