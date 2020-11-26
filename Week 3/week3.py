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
#       Update: alternatieve uitdrukking voor t, is die beter?
# - waarschijnlijk zijn er efficientere methodes voor FTotalOnAtoms en dat ding printen
# - wat Ruben zei in Teams over volgorde van atomen (in group X donderdag 19 nov)
#       maar dat komt dus volgende week

# BOND
def Vbond(r, k, r0):
    return(1/2*k*(r-r0)**2)

def Fbond(r, k, r0):
    return(-k*(r-r0))

def FBondOnAtoms(a, b, k, r0):
    r = np.linalg.norm(a-b)
    Fa = Fbond(r, k, r0)*(a-b)/np.linalg.norm(a-b)
    return(np.asarray([Fa, -Fa]))

# ANGLE
def Vangle(t, kt, t0):
    return 1/2*kt*(t-t0)**2

def Fangle(t, kt, t0):
    return -kt*(t-t0)

def FAngleOnAtoms(a, b, c, kt, t0):
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
    Fa = Fangle(t, kt, t0)/np.linalg.norm(a-b) * normalVecA/np.linalg.norm(normalVecA)
    Fc = Fangle(t, kt, t0)/np.linalg.norm(c-b) * normalVecC/np.linalg.norm(normalVecC)
    return(np.asarray([Fa, -Fa-Fc, Fc]))

def FTotalOnAtoms(a,b,c, k, r0, kt, t0):
    """Compute total forces on 3-body atom."""
    FBondAB = FBondOnAtoms(a,b,k,r0)
    FBondBC = FBondOnAtoms(b,c,k,r0)
    FAngle = FAngleOnAtoms(a,b,c,kt,t0)
    Fa = FBondAB[0] + FAngle[0]
    Fc = FBondBC[1] + FAngle[2]
    Fb = FBondAB[1] + FBondBC[0] + FAngle[1]
    return(np.asarray([Fb, Fa, Fc]))

# hydrogen example
types, xyzs = readXYZfile("HydrogenSingle.xyz", 0)
k = 24531/(10**2) # in kJ / (mol A^2)
r0 = 0.74 # in Angstrom

# print(FBondOnAtoms(xyzs[0],xyzs[1], k, r0))

# water example
types, xyzs = readXYZfile("WaterSingle.xyz", 0)
k = 502416/(10**2) # in kJ / (mol A^2)
r0 = 0.9572 # in Angstrom
kt = 628.02
t0 = np.deg2rad(104.52)

# print(FTotalOnAtoms(xyzs[1] , xyzs[0], xyzs[2], k, r0, kt, t0))



### WEEK 3 ###
# TODO: 
#   - Do integrators need a? Don't think so, if force and mass are known
#   - Verlet not yet working. Using two timesteps back seems to be a problem

def integratorEuler(x, v, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for this integrator. """ 
    # print(dt*FBondOnAtoms(x[0], x[1], k, r0)/m)
    """ Implementation of a single step for Euler integrator. """ 
    x = x + dt*v + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m
    v = v + dt*FBondOnAtoms(x[0], x[1], k, r0)/m
    return(x, v, a)


def integratorVerlet(x, x1, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for Verlet integrator. """ 
    x = 2*x - x1  + (dt**2) * FBondOnAtoms(x[0], x[1], k, r0)/m
    v = 1/(2*dt) * (x - x1)
    return(x, v, a)

def integratorVerlocity(x, x1, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for Verlet integrator. """ 
    x = 2*x - x1  + (dt**2) * FBondOnAtoms(x[0], x[1], k, r0)/m
    v = 1/(2*dt) * (x - x1)
    return(x, v, a)

# Add RK4?

# Generate a random velocity:
# first get a random unit vector (direction) 


# H2 example
time = 0
endTime = 1
types, x = readXYZfile("HydrogenSingle.xyz", 0)
k = 24531/(10**2) # in kJ / (mol A^2)
r0 = 0.74 # in Angstrom
u1 = np.random.uniform(size=3)
u2 = np.random.uniform(size=3)
u1 /= np.linalg.norm(u1) # normalize
u2 /= np.linalg.norm(u2)
v1 = 0.01*u1
v2 = 0.01*u2

# To make initial velocity zero, uncomment these lines:
# v1 = np.array([0,0,0])
# v2 = np.array([0,0,0])

v = np.asarray([v1,v2])
m = 1.00784
a = FBondOnAtoms(x[0],x[1], k, r0)/m
kt = 0
t0 = 0
dt = 0.1 * 2*np.pi*np.sqrt(m/k)
# Might need other estimate for larger molecules

with open("output.txt", "w") as outputFile: # clear output file 
        outputFile.write("")

# while(time<=endTime):
#     # print(x)
#     x, v, a = integratorEuler(x, v, a, m, k, r0, kt, t0, dt) 
#     time += dt
    
with open("output.xyz", "a") as outputFile:
    outputFile.write(f"{len(types)}\n")
    outputFile.write(f"This is a comment and the time is {time:5.4f}\n")
    for i, atom in enumerate(x):
        outputFile.write(f"{types[i]} {x[i,0]:10.5f} {x[i,1]:10.5f} {x[i,2]:10.5f}\n")

# Euler example: 
# while(time<=endTime):
#     print(x)
#     x, v, a = integratorEuler(x, v, a, m, k, r0, kt, t0, dt) 
#     time += dt
# print(x)

#Verlet example: (should use Euler as first step)
nrTimeSteps = int(endTime/dt)
types, xmin1 = readXYZfile("HydrogenSingle.xyz", 0)
x, v, a = integratorEuler(xmin1, v, a, m, k, r0, kt, t0, dt) 
for i in range(nrTimeSteps):
    print(x)
    x1_temp = xmin1 
    xmin1 = x # to store x[(i+1)-1] for next iteration
    x, v, a = integratorVerlet(x, x1_temp, a, m, k, r0, kt, t0, dt) 



































