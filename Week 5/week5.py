# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)

"""

import numpy as np

### WEEK 1 ###
# To do (misschien): geheugen-efficiënter maken
def readXYZnrAtoms(fileName): 
    """Reads number of atoms given in xyz file """
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
    
    atomTypes = firstColumn[(2+(2+nrOfAtoms)*timeStep):((2+nrOfAtoms)*(timeStep+1))]
    atomPositions = lines[(2+(2+nrOfAtoms)*timeStep):((2+nrOfAtoms)*(timeStep+1))]
        
    atomPositions = np.asarray(atomPositions).astype(np.float)
    return(atomTypes,atomPositions)

def distAtoms(positions):
    """ Computes distances between all atoms """
    diff = positions - positions[:,np.newaxis]
    dist = np.linalg.norm(diff,axis = 2)
    return(dist)


### WEEK 2 ###

# BOND
def Vbond(r, k, r0):
    """ Calculates harmonic bond potential """
    return(1/2*k*(r-r0)**2)

def Fbond(r, k, r0):
    """ Calculates bond force magnitude """
    return(-k*(r-r0))


# ANGLE
def Vangle(t, kt, t0):
    """ Calculates harmonic angular potential """
    return 1/2*kt*(t-t0)**2

def Fangle(t, kt, t0):
    """ Calculates angular force magnitude """
    return -kt*(t-t0)



### WEEK 3 updated ###
def integratorEulerNew(x, v, a):
    """ Implementation of a single step for Euler integrator. """ 
    x = x + dt*v + (dt**2)/2*a
    v = v + dt*a
    return(x, v, a)

def integratorVerletNew(x, x1, v, a): # not yet working since x1 is needed
    """ Implementation of a single step for Verlet integrator. """ 
    x = 2*x - x1  + (dt**2) * a
    v = 1/(2*dt) * (x - x1)
    return(x, v, a)

def integratorVerlocityNew(x, v, a):
    """ Implementation of a single step for Velocty Verlet integrator. """ 
    x_new = x + v*dt + (dt**2)/2*a
    a_new = computeForces(x_new, bonds, bondConstants, angles, angleConstants)/m[:,np.newaxis]
    v = v + dt/2*(a_new +a)
    return(x_new, v, a)

def integratorRK4New(x, v, a):
    """ Implementation of a single step for Runge-Kutta order 4 integrator. """ 
    x1 = x + dt*v + (dt**2)/2*a 
    v1 = dt*computeForces(x1, bonds, bondConstants, angles, angleConstants)/m[:,np.newaxis]
    x2 = x + dt/2*(v+v1/2) + (dt**2)/2*a 
    v2 = dt*computeForces(x2, bonds, bondConstants, angles, angleConstants)/m[:,np.newaxis]
    x3 = x + dt/2*(v+v2/2) + (dt**2)/2*a
    v3 = dt*computeForces(x3, bonds, bondConstants, angles, angleConstants)/m[:,np.newaxis]
    x4 = x + dt*(v+v3) + (dt**2)/2*a 
    v4 = dt*computeForces(x4, bonds, bondConstants, angles, angleConstants)/m[:,np.newaxis]
    
    v = v + (v1+2*v2+2*v3+v4)/6
    return(x1, v, a)

### WEEK 4 ###
# TODO: 
# - Fix integrators to work on only the input: x,v,a
# - Only need: readXYZfile, Fbond, Fangle and integrators from previous weeks

def readTopologyFile(fileNameTopology): 
    """Read a topology file."""
    with open(fileNameTopology, "r") as inputFile:
        lines = inputFile.readlines()
        
        ### ik dacht dat het handig was om aan te geven welke atomen in een molecuul zitten
        ### maar dat is tot nu to nog niet echt nodig gebleken
        nrOfMolecules = int(lines[0].split()[1])
        # print(nrOfMolecules)
    
        molecules = []
        for i in range(1,nrOfMolecules+1):
            molecules.append(list(map(int,lines[i].split())))
        # molecules = np.asarray(molecules) # no np array for molecules since length differs per entry
        # print(molecules)
        ###

        nrOfBonds = int(lines[nrOfMolecules+1].split()[1])
        bonds = []
        bondConstants = []
        for i in range(nrOfMolecules+2,nrOfMolecules+nrOfBonds+2):
            bonds.append([int(lines[i].split()[0]),int(lines[i].split()[1])])
            bondConstants.append([float(lines[i].split()[2]),float(lines[i].split()[3])])
        bonds = np.asarray(bonds).reshape(nrOfBonds,2)
        bondConstants = np.asarray(bondConstants).reshape(nrOfBonds,2)
    
        nrOfAngles = int(lines[nrOfMolecules+nrOfBonds+2].split()[1])
        angles = []
        angleConstants = []
        for i in range(nrOfMolecules+nrOfBonds+3,len(lines)):
            angles.append([int(lines[i].split()[0]),int(lines[i].split()[1]),int(lines[i].split()[2])])
            angleConstants.append([float(lines[i].split()[3]),float(lines[i].split()[4])])
        angles = np.asarray(angles).reshape(nrOfAngles,3)
        angleConstants = np.asarray(angleConstants).reshape(nrOfAngles,2)
        
        return(molecules, bonds, bondConstants, angles, angleConstants)


def computeForces(x, bonds, bondConstants, angles, angleConstants):
    """Caltulate forces in one go with help of topology file."""
    forces = np.zeros((len(types),3), dtype = float)
    # bonds
    if bonds.size > 0:
        r = np.linalg.norm(x[bonds[:,0]] - x[bonds[:,1]], axis = 1)
        Fbonds = Fbond(r, bondConstants[:,0], bondConstants[:,1])[:,np.newaxis]*(x[bonds[:,0]]-x[bonds[:,1]])/r[:,np.newaxis]
        np.add.at(forces, bonds[:,0], Fbonds)
        np.add.at(forces, bonds[:,1], -Fbonds)
    
    # angles 
    if angles.size > 0:
        atomLeft = x[angles[:,0]]
        atomMiddle = x[angles[:,1]]
        atomRight = x[angles[:,2]]
        dif1 = atomLeft-atomMiddle
        dif2 = atomRight-atomMiddle
        t = np.arccos(np.sum(dif1*dif2, axis = 1)/(np.linalg.norm(atomLeft-atomMiddle, axis = 1)*np.linalg.norm(atomRight-atomMiddle, axis = 1)))
        Fangles = Fangle(t, angleConstants[:,0], angleConstants[:,1])

        normalVec1 = np.cross(atomLeft-atomMiddle,np.cross(atomLeft-atomMiddle,atomRight-atomMiddle))
        normalVec2 = np.cross(atomMiddle-atomRight,np.cross(atomLeft-atomMiddle,atomRight-atomMiddle))

        FangleAtomLeft = Fangles[:,np.newaxis]/np.linalg.norm(atomLeft-atomMiddle, axis = 1)[:,np.newaxis] * normalVec1/np.linalg.norm(normalVec1, axis = 1)[:,np.newaxis]
        FangleAtomRight = Fangles[:,np.newaxis]/np.linalg.norm(atomRight-atomMiddle, axis = 1)[:,np.newaxis] * normalVec2/np.linalg.norm(normalVec2, axis = 1)[:,np.newaxis]
        FangleAtomMiddle = -FangleAtomLeft - FangleAtomRight
        
        np.add.at(forces, angles[:,0], FangleAtomLeft)
        np.add.at(forces, angles[:,1], FangleAtomMiddle)
        np.add.at(forces, angles[:,2], FangleAtomRight)
    return(forces)

# example
# two water and one hydrogen molecules
types, x = readXYZfile("MixedMolecules.xyz", 0)
molecules, bonds, bondConstants, angles, angleConstants = readTopologyFile("MixedMoleculesTopology.txt")

time_loc = 0
endTime = 2
dt = 0.001


massesDict = {'H': 1.00784, 'O': 15.9994, 'C': 12.0110}
m = np.vectorize(massesDict.get)(types)
# kan ook direct bij xyz-lezen


x_loc = x
u = np.random.uniform(size=3*len(types)).reshape((len(types),3)) # random starting velocity vector
u = u/np.linalg.norm(u,axis = 1)[:,np.newaxis] # normalize
v_loc = 0.1*u

with open("MixedMoleculesOutput.xyz", "w") as outputFile: 
    outputFile.write("") 
with open("MixedMoleculesOutput.xyz", "a") as outputFile:
    while (time_loc <= endTime) : 
        outputFile.write(f"{len(types)}\n")
        outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
        for i, atom in enumerate(x_loc):
            outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")  
            
        forces = computeForces(x_loc, bonds, bondConstants, angles, angleConstants)
        accel = forces / m[:,np.newaxis]
        x_loc, v_loc, a_loc = integratorEulerNew(x_loc, v_loc, accel)
        time_loc += dt
          









