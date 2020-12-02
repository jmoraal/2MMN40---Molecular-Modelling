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

# TODO
# - waarschijnlijk zijn er efficientere methodes voor FTotalOnAtoms en dat ding printen
# - wat Ruben zei in Teams over volgorde van atomen (in group X donderdag 19 nov)
#       maar dat komt dus volgende week

# BOND
def Vbond(r, k, r0):
    """ Calculates harmonic bond potential """
    return(1/2*k*(r-r0)**2)

def Fbond(r, k, r0):
    """ Calculates bond force magnitude """
    return(-k*(r-r0))

def FBondOnAtoms(a, b, k, r0):
    """ Calculates bond force (with direction) """
    r = np.linalg.norm(a-b)
    Fa = Fbond(r, k, r0)*(a-b)/np.linalg.norm(a-b)
    return(np.asarray([Fa, -Fa]))

# ANGLE
def Vangle(t, kt, t0):
    """ Calculates harmonic angular potential """
    return 1/2*kt*(t-t0)**2

def Fangle(t, kt, t0):
    """ Calculates angular force magnitude """
    return -kt*(t-t0)


def FAngleOnAtoms(o, h1, h2, kt, t0):
    """ Compute angular forces on 3-body atom.
    
    Mind the order of the arguments. The first argument is supposed to be the middle atom (O in water).
    INPUT: positions of atoms a,b,c
    OUTPUT: angular force acting on each of the atoms
    """
    t = np.arccos(np.dot((h1-o),(h2-o))/(np.linalg.norm(h1-o)*np.linalg.norm(h2-o)))
    """ Alternative computation of t using cosine rule:
    ab = np.linalg.norm(a-b)
    bc = np.linalg.norm(c-b)
    ac = np.linalg.norm(a-c)
    t = np.arccos((ab**2 + bc**2 - ac**2)/(2*ab*bc))
    """
    normalVech1 = np.cross(h1-o,np.cross(h1-o,h2-o))
    normalVech2 = np.cross(o-h2,np.cross(h1-o,h2-o))
    Fh1 = Fangle(t, kt, t0)/np.linalg.norm(h1-o) * normalVech1/np.linalg.norm(normalVech1)
    Fh2 = Fangle(t, kt, t0)/np.linalg.norm(h2-o) * normalVech2/np.linalg.norm(normalVech2)
    return(np.asarray([-Fh1-Fh2, Fh1, Fh2]))

def FTotalOnAtoms(centre, outer1, outer2, k, r0, kt, t0):
    """ Compute total forces on 3-body atom. """
    FBondouter1centre = FBondOnAtoms(outer1,centre,k,r0)
    FBondcentreouter2 = FBondOnAtoms(centre,outer2,k,r0)
    FAngle = FAngleOnAtoms(centre,outer1,outer2,kt,t0)
    Fouter1 = FBondouter1centre[0] + FAngle[1]
    Fouter2 = FBondcentreouter2[1] + FAngle[2]
    Fcentre = FBondouter1centre[1] + FBondcentreouter2[0] + FAngle[0]
    return(np.asarray([Fcentre, Fouter1, Fouter2]))

# hydrogen example
def hydrogenForcesExample():
    types, xyzs = readXYZfile("HydrogenSingle.xyz", 0)
    k = 24531/(10**2) # in kJ / (mol A^2)
    r0 = 0.74 # in Angstrom
    print(FBondOnAtoms(xyzs[0],xyzs[1], k, r0))

# water example
def waterForcesExample():
    types, xyzs = readXYZfile("WaterSingle.xyz", 0)
    k = 502416/(10**2) # in kJ / (mol A^2)
    r0 = 0.9572 # in Angstrom
    kt = 628.02
    t0 = np.deg2rad(104.52)
    print(FTotalOnAtoms(xyzs[0], xyzs[1], xyzs[2], k, r0, kt, t0))

# waterForcesExample()

### WEEK 3 ###
# TODO: 
#   - Do integrators need a? Don't think so at the moment
#   - Generalise integrators (remove case distinction 2/3-body)
#   - Nicer implementation for 'x_loc = x' construction?

def integratorEuler(x, v, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for Euler integrator. """ 
    if len(types) == 2:
        x = x + dt*v + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m
        v = v + dt*FBondOnAtoms(x[0], x[1], k, r0)/m
    elif len(types) == 3:
        x = x + dt*v + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m
        v = v + dt*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m
    return(x, v, a)

def integratorEulerNew(x, v, a):
    """ Implementation of a single step for Euler integrator. """ 
    x = x + dt*v + (dt**2)/2*a
    v = v + dt*a
    return(x, v, a)

# Verlet was also implemented out of curiosity, but turned out to be inconvenient to handle. 
def integratorVerlet(x, x1, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for Verlet integrator. """ 
    if len(types) == 2:
        x = 2*x - x1  + (dt**2) * FBondOnAtoms(x[0], x[1], k, r0)/m
        v = 1/(2*dt) * (x - x1)
    elif len(types) == 3:
        x = 2*x - x1  + (dt**2) * FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m
        v = 1/(2*dt) * (x - x1)
    return(x, v, a)


def integratorVerlocity(x, v, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for Velocty Verlet integrator. """ 
    if len(types) == 2:
        x_new = x + v*dt + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m 
        v = v + dt/2*(FBondOnAtoms(x[0], x[1], k, r0)/m + FBondOnAtoms(x_new[0], x_new[1], k, r0)/m)
    elif len(types) == 3:
        x_new = x + v*dt + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m 
        v = v + dt/2*(FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m + FTotalOnAtoms(x_new[0], x_new[1], x_new[2], k, r0, kt, t0)/m)
    
    return(x_new, v, a)

def integratorVerlocityNew(x, v, a):
    """ Implementation of a single step for Velocty Verlet integrator. """ 
    x_new = x + v*dt + (dt**2)/2*a
    a_new = 1
    v = v + dt/2*(a_new +a)
    return(x_new, v, a)

def integratorRK4(x, v, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for Runge-Kutta order 4 integrator. """ 
    if len(types) == 2:
        x1 = x + dt*v + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m 
        v1 = dt*FBondOnAtoms(x1[0], x1[1], k, r0)/m 
        x2 = x + dt/2*(v+v1/2) + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m 
        v2 = dt*FBondOnAtoms(x2[0], x2[1], k, r0)/m 
        x3 = x + dt/2*(v+v2/2) + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m 
        v3 = dt*FBondOnAtoms(x3[0], x3[1], k, r0)/m 
        x4 = x + dt*(v+v3) + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m 
        v4 = dt*FBondOnAtoms(x4[0], x4[1], k, r0)/m 
        
        x = x + dt*v + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m 
        v = v + (v1+2*v2+2*v3+v4)/6
    elif len(types) == 3:
        x1 = x + dt*v + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m 
        v1 = dt*FTotalOnAtoms(x1[0], x1[1], x1[2], k, r0, kt, t0)/m 
        x2 = x + dt/2*(v+v1/2) + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m 
        v2 = dt*FTotalOnAtoms(x2[0], x2[1], x2[2], k, r0, kt, t0)/m 
        x3 = x + dt/2*(v+v2/2) + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m 
        v3 = dt*FTotalOnAtoms(x3[0], x3[1], x3[2], k, r0, kt, t0)/m 
        x4 = x + dt*(v+v3) + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m 
        v4 = dt*FTotalOnAtoms(x4[0], x4[1], x4[2], k, r0, kt, t0)/m 
        
        x = x + dt*v + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m 
        v = v + (v1+2*v2+2*v3+v4)/6
    return(x, v, a)

def integratorRK4New(x, v, a, a1, a2, a3, a4):
    """ Implementation of a single step for Runge-Kutta order 4 integrator. """ 
    x1 = x + dt*v + (dt**2)/2*a 
    v1 = dt*a1 
    x2 = x + dt/2*(v+v1/2) + (dt**2)/2*a 
    v2 = dt*a2
    x3 = x + dt/2*(v+v2/2) + (dt**2)/2*a
    v3 = dt*a3
    x4 = x + dt*(v+v3) + (dt**2)/2*a 
    v4 = dt*a4 
    
    x = x + dt*v + (dt**2)/2*a
    v = v + (v1+2*v2+2*v3+v4)/6
    return(x, v, a)

def setParametersH2 (velocityZero =False) :
    """ Sets example parameters for a single hydrogen molecule """
    global time, endTime, types, x, k, r0, v1, v2, v, m, a, kt, t0, dt
 
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
    
    if velocityZero : 
        v1 = np.array([0,0,0])
        v2 = np.array([0,0,0])
    
    v = np.asarray([v1,v2])
    m = 1.00784
    a = FBondOnAtoms(x[0],x[1], k, r0)/m
    kt = 0
    t0 = 0
    dt = 0.1 * 2*np.pi*np.sqrt(m/k)
    # Might need other estimate for larger molecules


def setParametersH2O (velocityZero =False) :
    """ Sets example parameters for a single water molecule """
    global time, endTime, types, x, k, r0, v1, v2, v, m, a, kt, t0, dt
 
    time = 0
    endTime = 3
    types, x = readXYZfile("WaterSingle.xyz", 0)
    k = 502416/(10**2) # in kJ / (mol A^2)
    r0 = 0.9572 # in Angstrom
    kt = 628.02
    t0 = np.deg2rad(104.52)
    u1 = np.random.uniform(size=3)
    u2 = np.random.uniform(size=3)
    u3 = np.random.uniform(size=3)
    u1 /= np.linalg.norm(u1) # normalize
    u2 /= np.linalg.norm(u2)
    u3 /= np.linalg.norm(u3)
    v1 = 0.01*u1 
    v2 = 0.01*u2
    v3 = 0.01*u3
    # effect only significant for 0.1*u not for 0.01*u
    
    if velocityZero : 
        v1 = np.array([0,0,0])
        v2 = np.array([0,0,0])
        v3 = np.array([0,0,0])
    
    v = np.asarray([v1,v2,v3])
    m = np.asarray([15.999,1.00784,1.00784])
    a = FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m
    dt = 0.1 * 2*np.pi*np.sqrt(np.amin(m)/k)
    # Might need other estimate for larger molecules?


def writeExampleToXYZ(molecule, integrator, filename, velocityZero =False):
    """ Write examples of molecule movements to xyz
    
    Important: does not work for (standard) Verlet! Only Euler, Verlocity, RK4
    INPUT: molecule type (string), integrator (function), and possibly whether initial velocity should be zero
    OUTPUT: movement of given molecule as computed by specified integrator, as xyz-fle
    """
    if (molecule == 'H2') : 
        setParametersH2(velocityZero)   
    elif (molecule == 'H2O'):
        setParametersH2O(velocityZero) 
    else :
        print("No settings for this molecule")
        # Maybe a nicer error catch?
        
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time

    with open(filename, "w") as outputFile: # clear output file 
            outputFile.write("") #Can we let the file name depend on molecule and integrator?
    
    
    with open(filename, "a") as outputFile:
        while (time_loc <= endTime) : 
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")  
            x_loc, v_loc, a_loc = integrator(x_loc, v_loc, a_loc, m, k, r0, kt, t0, dt) 
            time_loc += dt

    #return x_loc, v_loc, a_loc
 
# writeExampleToXYZ('H2', integratorEuler, "EulerH2Example.xyz")
# writeExampleToXYZ('H2', integratorVerlocity, "VerlocityH2Example.xyz")
# writeExampleToXYZ('H2', integratorRK4, "RK4H2Example.xyz")
# writeExampleToXYZ('H2O', integratorEuler, "EulerH2OExample.xyz")
# writeExampleToXYZ('H2O', integratorVerlocity, "VerlocityH2OExample.xyz")
# writeExampleToXYZ('H2O', integratorRK4, "RK4H2OExample.xyz")

def eulerNewExample(molecule, filename, velocityZero =False):
    """ Not working yet for H2O """
    if (molecule == 'H2') : 
        setParametersH2(velocityZero)   
    elif (molecule == 'H2O'):
        setParametersH2O(velocityZero) 
    else :
        print("No settings for this molecule")
        # Maybe a nicer error catch?
        
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time

    with open(filename, "w") as outputFile: # clear output file 
            outputFile.write("") #Can we let the file name depend on molecule and integrator?
    
    
    with open(filename, "a") as outputFile:
        while (time_loc <= endTime) : 
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")  
            
            if len(types) == 2:
                forces = FBondOnAtoms(x_loc[0], x_loc[1], k, r0)
            elif len(types) == 3:
                forces = FTotalOnAtoms(x_loc[0], x_loc[1], x_loc[2], k, r0, kt, t0)
            
            accel = forces / m
            #accel = np.true_divide(forces, m[:,np.newaxis])
            
            x_loc, v_loc, a_loc = integratorEulerNew(x_loc, v_loc, accel)
            time_loc += dt

    #return x_loc, v_loc, a_loc
eulerNewExample('H2', "EulerH2OExample.xyz")

# def verlocityNewExample(molecule, filename, velocityZero =False):
#     """ Write examples of molecule movements to xyz
    
#     NOT WORKING YET
#     """
#     if (molecule == 'H2') : 
#         setParametersH2(velocityZero)   
#     elif (molecule == 'H2O'):
#         setParametersH2O(velocityZero) 
#     else :
#         print("No settings for this molecule")
#         # Maybe a nicer error catch?
        
#     x_loc = x
#     v_loc = v
#     a_loc = a
#     time_loc = time

#     with open(filename, "w") as outputFile: # clear output file 
#             outputFile.write("") #Can we let the file name depend on molecule and integrator?
    
    
#     with open(filename, "a") as outputFile:
#         while (time_loc <= endTime) : 
#             outputFile.write(f"{len(types)}\n")
#             outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
#             for i, atom in enumerate(x_loc):
#                 outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")  
            
#             x_loc, v_loc, a_loc = integrator(x_loc, v_loc, a_loc, m, k, r0, kt, t0, dt) 
#             time_loc += dt

#     #return x_loc, v_loc, a_loc



### Trying to read topology and work with different types of molecules simultaneously ###
types, x = readXYZfile("MixedMolecules.xyz", 0)
with open("MixedMoleculesTopology.txt", "r") as inputFile:
    lines = inputFile.readlines()
    nrOfMolecules = int(lines[0].split()[1])
    print(nrOfMolecules)
    
    molecules = []
    for i in range(1,nrOfMolecules+1):
        molecules.append(list(map(int,lines[i].split())))
    # molecules = np.asarray(molecules) # no np array for molecules since length differs per entry
    # print(molecules)

    nrOfBonds = int(lines[nrOfMolecules+1].split()[1])
    bonds = []
    bondConstants = []
    for i in range(nrOfMolecules+2,nrOfMolecules+nrOfBonds+2):
        bonds.append([int(lines[i].split()[0]),int(lines[i].split()[1])])
        bondConstants.append([float(lines[i].split()[2]),float(lines[i].split()[3])])
    bonds = np.asarray(bonds).reshape(nrOfBonds,2)
    bondConstants = np.asarray(bondConstants).reshape(nrOfBonds,2)
    # print(bonds)
    # print(bondConstants)
    
    nrOfAngles = int(lines[nrOfMolecules+nrOfBonds+2].split()[1])
    angles = []
    angleConstants = []
    for i in range(nrOfMolecules+nrOfBonds+3,len(lines)):
        print(i)
        angles.append([int(lines[i].split()[0]),int(lines[i].split()[1]),int(lines[i].split()[2])])
        angleConstants.append([float(lines[i].split()[3]),float(lines[i].split()[4])])
    angles = np.asarray(angles).reshape(nrOfAngles,3)
    angleConstants = np.asarray(angleConstants).reshape(nrOfAngles,2)
    # print(angles)
    # print(angleConstants)
    
with open("MixedMoleculesOutput.xyz", "w") as outputFile: 
    outputFile.write("") 
# with open("WaterDoubleOutput", "a") as outputFile:
    # while (time_loc <= endTime) : 
    #     outputFile.write(f"{len(types)}\n")
    #     outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
          














