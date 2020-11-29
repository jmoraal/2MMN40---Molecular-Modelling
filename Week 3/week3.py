# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981

NOTE: using tab for indentation
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

# def FAngleOnAtoms(a, b, c, kt, t0):
#     """Compute angular forces on 3-body atom.
    
#     Mind the order of the arguments. The middle argument is supposed to be the middle atom (O in water).
#     INPUT: positions of atoms a,b,c
#     OUTPUT: angular force acting on each of the atoms
#     """
#     t = np.arccos(np.dot((a-b),(c-b))/(np.linalg.norm(a-b)*np.linalg.norm(c-b)))
#     """ Alternative computation of t using cosine rule:
#     ab = np.linalg.norm(a-b)
#     bc = np.linalg.norm(c-b)
#     ac = np.linalg.norm(a-c)
#     t = np.arccos((ab**2 + bc**2 - ac**2)/(2*ab*bc))
#     """
#     normalVecA = np.cross(a-b,np.cross(a-b,c-b))
#     print("normalVecA")
#     print(normalVecA)
    
#     normalVecC = np.cross(b-c,np.cross(a-b,c-b))
#     print("normalVecC")
#     print(normalVecC)
#     Fa = Fangle(t, kt, t0)/np.linalg.norm(a-b) * normalVecA/np.linalg.norm(normalVecA)
#     Fc = Fangle(t, kt, t0)/np.linalg.norm(c-b) * normalVecC/np.linalg.norm(normalVecC)
#     return(np.asarray([Fa, -Fa-Fc, Fc]))

# def FTotalOnAtoms(a,b,c, k, r0, kt, t0):
#     """Compute total forces on 3-body atom."""
#     FBondAB = FBondOnAtoms(a,b,k,r0)
#     FBondBC = FBondOnAtoms(b,c,k,r0)
#     FAngle = FAngleOnAtoms(a,b,c,kt,t0)
#     print(FAngle)
#     Fa = FBondAB[0] + FAngle[0]
#     Fc = FBondBC[1] + FAngle[2]
#     Fb = FBondAB[1] + FBondBC[0] + FAngle[1]
#     return(np.asarray([Fb, Fa, Fc]))

def FAngleOnAtoms(o, h1, h2, kt, t0):
    """Compute angular forces on 3-body atom.
    
    Mind the order of the arguments. The middle argument is supposed to be the middle atom (O in water).
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

def FTotalOnAtoms(o, h1, h2, k, r0, kt, t0):
    """Compute total forces on 3-body atom."""
    FBondh1o = FBondOnAtoms(h1,o,k,r0)
    FBondoh2 = FBondOnAtoms(o,h2,k,r0)
    FAngle = FAngleOnAtoms(o,h1,h2,kt,t0)
    Fh1 = FBondh1o[0] + FAngle[1]
    Fh2 = FBondoh2[1] + FAngle[2]
    Fo = FBondh1o[1] + FBondoh2[0] + FAngle[0]
    return(np.asarray([Fo, Fh1, Fh2]))

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
#   - Do integrators need a? Don't think so, if force and mass are known
#   - Verlet not yet working. Using two timesteps back seems to be a problem
#   - Add anglular forces into integrators

def integratorEuler(x, v, a, m, k, r0, kt, t0, dt):
    """ Implementation of a single step for this integrator. """ 
    if len(types) == 2:
        x = x + dt*v + (dt**2)/2*FBondOnAtoms(x[0], x[1], k, r0)/m
        v = v + dt*FBondOnAtoms(x[0], x[1], k, r0)/m
    elif len(types) == 3:
        x = x + dt*v + (dt**2)/2*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m
        v = v + dt*FTotalOnAtoms(x[0], x[1], x[2], k, r0, kt, t0)/m
    return(x, v, a)


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


# Add RK4?

# Generate a random velocity:
# first get a random unit vector (direction) 


# H2 example
def setParametersH2 (velocityZero =False) :
    global time, endTime
    global types, x
    global k, r0
    global v1, v2, v
    global m
    global a
    global kt
    global t0
    global dt  
 
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



def EulerH2Example(velocityZero =False) : 
    setParametersH2(velocityZero)
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time
    
    with open("EulerH2Example.xyz", "w") as outputFile: # clear output file 
        outputFile.write("")
    
    while(time_loc<=endTime):
        with open("EulerH2Example.xyz", "a") as outputFile:
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")
                
        x_loc, v_loc, a_loc = integratorEuler(x_loc, v_loc, a_loc, m, k, r0, kt, t0, dt) 
        time_loc += dt

    # print(x_loc)
    
    return x_loc, v_loc, a_loc


def VerletH2Example(velocityZero =False) : 
    setParametersH2(velocityZero)
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time
    
    types, xmin1 = readXYZfile("HydrogenSingle.xyz", 0)
    x_loc, v_loc, a_loc = integratorEuler(xmin1, v, a, m, k, r0, kt, t0, dt) 
    
    with open("VerletH2Example.xyz", "w") as outputFile: # clear output file 
        outputFile.write("")
    
    while(time_loc<=endTime):
        with open("VerletH2Example.xyz", "a") as outputFile:
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")
                
        # print(x_loc)
        x1_temp = xmin1 
        xmin1 = x_loc # to store x[(i+1)-1] for next iteration
        x_loc, v_loc, a_loc = integratorVerlet(x_loc, x1_temp, a_loc, m, k, r0, kt, t0, dt) 
        time_loc += dt
        
    return x_loc, v_loc, a_loc

def VerlocityH2Example(velocityZero =False) : 
    setParametersH2(velocityZero)
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time
    with open("VerlocityH2Example.xyz", "w") as outputFile: # clear output file 
        outputFile.write("")
    
    while(time_loc<=endTime):
        with open("VerlocityH2Example.xyz", "a") as outputFile:
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")
                
        # print(x_loc)
        x_loc, v_loc, a_loc = integratorVerlocity(x_loc, v_loc, a_loc, m, k, r0, kt, t0, dt)
        time_loc += dt
    
    return x_loc, v_loc, a_loc


EulerH2Example(velocityZero = False)    
VerletH2Example(velocityZero = False)
VerlocityH2Example(velocityZero = False)


# H2O example
def setParametersH2O (velocityZero =False) :
    global time, endTime
    global types, x
    global k, r0
    global v1, v2, v
    global m
    global a
    global kt
    global t0
    global dt  
 
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
    # Might need other estimate for larger molecules

def EulerH2OExample(velocityZero =False) : 
    setParametersH2O(velocityZero)
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time
    
    with open("EulerH2OExample.xyz", "w") as outputFile: # clear output file 
        outputFile.write("")
    
    while(time_loc<=endTime):
        with open("EulerH2OExample.xyz", "a") as outputFile:
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")
                
        x_loc, v_loc, a_loc = integratorEuler(x_loc, v_loc, a_loc, m, k, r0, kt, t0, dt) 
        time_loc += dt
        
    # return x_loc, v_loc, a_loc

def VerletH2OExample(velocityZero =False) : 
    setParametersH2O(velocityZero)
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time
    
    types, xmin1 = readXYZfile("WaterSingle.xyz", 0)
    x_loc, v_loc, a_loc = integratorEuler(xmin1, v, a, m, k, r0, kt, t0, dt) 
    
    with open("VerletH2OExample.xyz", "w") as outputFile: # clear output file 
        outputFile.write("")
    
    while(time_loc<=endTime):
        with open("VerletH2OExample.xyz", "a") as outputFile:
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")
                
        # print(x_loc)
        x1_temp = xmin1 
        xmin1 = x_loc # to store x[(i+1)-1] for next iteration
        x_loc, v_loc, a_loc = integratorVerlet(x_loc, x1_temp, a_loc, m, k, r0, kt, t0, dt) 
        time_loc += dt
    
    # return x_loc, v_loc, a_loc

def VerlocityH2OExample(velocityZero =False) : 
    setParametersH2O(velocityZero)
    x_loc = x
    v_loc = v
    a_loc = a
    time_loc = time
    with open("VerlocityH2OExample.xyz", "w") as outputFile: # clear output file 
        outputFile.write("")
    
    while(time_loc<=endTime):
        with open("VerlocityH2OExample.xyz", "a") as outputFile:
            outputFile.write(f"{len(types)}\n")
            outputFile.write(f"This is a comment and the time is {time_loc:5.4f}\n")
            for i, atom in enumerate(x_loc):
                outputFile.write(f"{types[i]} {x_loc[i,0]:10.5f} {x_loc[i,1]:10.5f} {x_loc[i,2]:10.5f}\n")
                
        # print(x_loc)
        x_loc, v_loc, a_loc = integratorVerlocity(x_loc, v_loc, a_loc, m, k, r0, kt, t0, dt)
        time_loc += dt
        
    return x_loc, v_loc, a_loc


EulerH2OExample(False)
VerletH2OExample(False)
VerlocityH2OExample(False)



















