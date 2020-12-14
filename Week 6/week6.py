# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)

"""

import numpy as np
#from itertools import chain

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
    massesDict = {'H': 1.00784, 'O': 15.9994, 'C': 12.0110}
    m = np.vectorize(massesDict.get)(atomTypes)
    return(atomTypes, atomPositions ,m)

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
#Note: can also be used for harmonic dihedral potential

def Fangle(t, kt, t0):
    """ Calculates angular force magnitude """
    return -kt*(t-t0)

def Vdihedral(t, k0, t0, C):
    """ Calculates periodic dihedral potential for given forcefield constants C """
    psi = t - np.pi
    return 0.5*(C[1]*(1+np.cos(psi)) + C[2]*(1-np.cos(2*psi)) + C[3]*(1+np.cos(3*psi)) + C[4]*(1-np.cos(4*psi)))

def Fdihedral(t, k0, t0, C):
    """ Calculates periodic dihedral force magnitude """
    psi = t - np.pi
    return 0.5*(-C[1]*np.sin(psi) + 2*C[2]*np.sin(2*psi) - 3*C[3]*np.sin(3*psi) + 4*C[4]*np.sin(4*psi))

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
        
        nrOfMolecules = int(lines[0].split()[1])
        # print(nrOfMolecules)
    
        molecules = []
        for i in range(1,nrOfMolecules+1):
            molecules.append(list(map(int,lines[i].split())))
        
        # create boolean matrix to indicate which atoms are part of same molecule:
        # Probably can still be done more efficiently?
        notInSameMolecule = np.ones((len(types), len(types)), dtype=bool)
        for i,mol in enumerate(molecules):
            for j,at in enumerate(mol):
                for k,at2 in enumerate(mol):
                    # print(at, at2)
                    notInSameMolecule[at,at2] = False
        
                
        
        # for i,atom in enumerate(types):
        #     for j,atom2 in enumerate(types):
        #         print(np.where(molecules == i)[0])
        #         if np.where(molecules == i)[0] == np.where(molecules == j)[0]:
        #             notInSameMolecule[i,j] = False
        
        # pad = len(maxa(molecules, key=len)) # add padding to create np array 
        # molecules = np.rray([i + [-1]*(pad-len(i)) for i in molecules])
        
        # another option for the notInSameMolecule is the nonAdjacencyList:
        # nonAdjacencyList = []
        # for i in range(0,len(molecules)):
        #     for j in range(0,len(molecules[i])):
        #         nonAdjacencyList.append(list(chain(*(molecules[:i] + molecules[i+1:]))))

        nrOfBonds = int(lines[nrOfMolecules+1].split()[1])
        bonds = []
        bondConstants = []
        for i in range(nrOfMolecules+2, nrOfMolecules+nrOfBonds+2):
            bonds.append([int(lines[i].split()[0]),int(lines[i].split()[1])])
            bondConstants.append([float(lines[i].split()[2]),float(lines[i].split()[3])])
        bonds = np.asarray(bonds).reshape(nrOfBonds,2)
        bondConstants = np.asarray(bondConstants).reshape(nrOfBonds,2)
    
        nrOfAngles = int(lines[nrOfMolecules+nrOfBonds+2].split()[1])
        angles = []
        angleConstants = []
        for i in range(nrOfMolecules+nrOfBonds+3, nrOfMolecules+nrOfBonds+nrOfAngles+3):
            angles.append([int(lines[i].split()[0]),int(lines[i].split()[1]),int(lines[i].split()[2])])
            angleConstants.append([float(lines[i].split()[3]),float(lines[i].split()[4])])
        angles = np.asarray(angles).reshape(nrOfAngles,3)
        angleConstants = np.asarray(angleConstants).reshape(nrOfAngles,2)
        
        sigma = []
        epsilon = []
        for i in range(nrOfMolecules+nrOfBonds+nrOfAngles+4, len(lines)):
            sigma.append(float(lines[i].split()[0]))
            epsilon.append(float(lines[i].split()[1]))
        sigma = np.asarray(sigma)
        epsilon = np.asarray(epsilon)
        
        return(notInSameMolecule, bonds, bondConstants, angles, angleConstants, sigma, epsilon)



def LennardJonesInter(sigma,eps,a,b,r):
    """computes LJ inter-molecular force
    
    Maybe calculate r inside of function?
    Not used at moment, but maybe can still be generalised to arrays instead of single calculations"""
    #r = np.linalg.norm()
    epsilon = np.sqrt(eps[a]*eps[b])
    sigma = 0.5*(sigma[a] + sigma[b])
    return 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

def distAtomsPBC(positions, boxSize):
    """ Computes distances between all atoms, with boundaries
    
    Not entirely sure this is correct!
    """
    # Old (but correct)
    # diff = abs(positions - positions[:,np.newaxis]) % boxSize
    # does not have right direction vector
    
    diff = positions - positions[:,np.newaxis] 
    diff = diff % (0.5*boxSize)
    # idea: if dist > 0.5*boxsize in some direction (x, y or z), then there is a closer copy. 
    # subtracting 0.5*boxsize in every direction where it is too large yields direction vector to closest neighbour
    # Still check correctness!
    dist = np.linalg.norm(diff,axis = 2)
    return(dist) #diff,
# Direction and distance are usually both needed, right? 
# Could also just return difference vector and do distance calculation elsewhere

def computeForces(x, bonds, bondConstants, angles, angleConstants, sigma, epsilon, boxSize = np.infty):
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
        
    dihedrals = [] # only there to prevent error warnings while implementing
    # dihedrals
    if dihedrals.size > 0:
        dif1 = x[dihedrals[:,0]] - x[dihedrals[:,1]]
        difCommon = x[dihedrals[:,1]] - x[dihedrals[:,2]]
        dif2 = x[dihedrals[:,2]] - x[dihedrals[:,3]]
        
        normalVec1 = np.cross(dif1,difCommon)
        normalVec2 = np.cross(-difCommon,dif2)
        
        theta = np.arccos(np.sum(normalVec1*normalVec2, axis = 1)/(np.linalg.norm(normalVec1, axis = 1)*np.linalg.norm(normalVec2, axis = 1)))
        
        psi = theta - np.pi
        
        #voorbereiding is gedaan, nu moeten (analoog aan angles) de krachten zelf nog berekend worden met Fdihedrals() en richting gegeven worden
    
    # Lennard Jones forces
    if sigma.size > 0:
        dist = distAtomsPBC(x,boxSize)  
        #U = np.zeros((len(types), 3))
    
        e = np.sqrt(epsilon*epsilon[:,np.newaxis])
        s = 0.5*(sigma+sigma[:,np.newaxis])
        
        frac = np.divide(s, dist, out=np.zeros_like(s), where=dist!=0) # avoid division by 0
        U = 4*e*(frac**12 - frac**6)
        L = x-x[:,np.newaxis]
        
        U = U*notInSameMolecule # these forces do not apply on atoms in the same molecule!
        
        U = np.repeat(U, 3).reshape(len(types), len(types), 3)
        
        forces += np.sum(U*L, axis = 1)
    
    # hieronder de LJ forces met forloops
    # dist = distAtomsPBC(x,boxSizd)    
    # if sigma.size > 0:
    #     f = np.zeros((len(types), 3))
    #     for i,atom in enumerate(types):
    #         for j,atom2 in enumerate(types):
    #             if np.where(molecules == i)[0] != np.where(molecules == j)[0]:
    #                 e = np.sqrt(epsilon[i]*epsilon[j])
    #                 s = 0.5*(sigma[i] + sigma[j])
    #                 r = dist[i,j]
    #                 U = 4*e*((s/r)**12 - (s/r)**6)
    #                 # r = dist[i,j]
    #                 # U = LennardJonesInter(sigma, epsilon, i, j, r)
    #                 if U != 0:
    #                     f[i] += U*(x[j] - x[i])                 
        
    #     forces = forces + f 
    return(forces)

# example
# two water and one hydrogen molecules
types, x, m = readXYZfile("MixedMolecules.xyz", 0)
notInSameMolecule, bonds, bondConstants, angles, angleConstants, sigma, epsilon = readTopologyFile("MixedMoleculesTopology.txt")

time_loc = 0
endTime = 0.00
dt = 0.001

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
            
        forces = computeForces(x_loc, bonds, bondConstants, angles, angleConstants, sigma, epsilon)
        accel = forces / m[:,np.newaxis]
        x_loc, v_loc, a_loc = integratorEulerNew(x_loc, v_loc, accel)
        time_loc += dt
          



### WEEK 5 ###
# TODO: 
# - Implement Lennard-Jones Potential & forces
# - Add periodic boundary conditions


#Neighbourlists: 

# Adjacency matrix: sparse, so memory-inefficient
def neighbourMatrix(positions,cutoff): 
    """ returns adjacancy matrix for atoms closer than cutoff """
    return (distAtoms(positions) < cutoff)


def neighbourList(positions,cutoff): 
    """ returns list of pairs of atoms with distance < cutoff
    
    Note: still contains duplicates and pairs (x,x)
    Also, there is probably a more direct way?"""
    neighMatrix = neighbourMatrix(positions,cutoff)
    neighPairs = np.matrix.transpose(np.array(np.where(neighMatrix)))
    return neighPairs


#atom-wise projection into box: 
def coordProjectToBox(x, boxSize):
    """Projects all coordinates into cubic box of size boxSize"""
    return (x % boxSize)
#TODO: generalise to non-cubic shapes? 
# if not, then the function above is not necessary



# example (disclaimer: not yet fully correct, just to check methods so far are working)
# two water and one hydrogen molecules
types, x, m = readXYZfile("MixedMolecules.xyz", 0)
atomsInOtherMolecules, bonds, bondConstants, angles, angleConstants, sigma, epsilon = readTopologyFile("MixedMoleculesTopology.txt")

boxSizeExample = 10 # TODO: yet to choose meaningful value
cutoff = 0.5*boxSizeExample

time = 0
endTime = 0.001
dt = 0.001

u = np.random.uniform(size=3*len(types)).reshape((len(types),3)) # random starting velocity vector
u = u/np.linalg.norm(u,axis = 1)[:,np.newaxis] # normalize
v = 0.1*u

with open("MixedMoleculesPBCOutput.xyz", "w") as outputFile: 
    outputFile.write("") 
with open("MixedMoleculesPBCOutput.xyz", "a") as outputFile:
    while (time <= endTime) : 
        outputFile.write(f"{len(types)}\n")
        outputFile.write(f"This is a comment and the time is {time:5.4f}\n")
        for i, atom in enumerate(x):
            outputFile.write(f"{types[i]} {x[i,0]:10.5f} {x[i,1]:10.5f} {x[i,2]:10.5f}\n")  
            
        neighList = neighbourList(x, 0.5*boxSizeExample)
        #yet to make sure LJ is computed correctly for neighbours from different boxes
        x = x % boxSizeExample
        forces = computeForces(x, bonds, bondConstants, angles, angleConstants, sigma, epsilon, boxSize = boxSizeExample)
        accel = forces / m[:,np.newaxis]
        x, v, a = integratorEulerNew(x, v, accel)
        time += dt