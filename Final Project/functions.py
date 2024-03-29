# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 21:19:51 2021

@author: s161981
"""

import numpy as np


### READ INPUT ###

# XYZ:
def readXYZnrAtoms(fileName): 
    """Reads number of atoms given in xyz file """
    with open(fileName, "r") as inputFile:
        firstString = inputFile.readline().split()
        nrOfAtoms = int(firstString[0])
    return nrOfAtoms

# TODO (misschien): geheugen-efficiënter maken
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
    massesDict = {'H': 1.0080, 'O': 15.9994, 'C': 12.0110} # in Dalton; #TODO: better at another place
    m = np.vectorize(massesDict.get)(atomTypes)
    return(atomTypes, atomPositions ,m)

# TOPOLOGY
def readTopologyFile(fileNameTopology): 
    """Read a topology file."""
    with open(fileNameTopology, "r") as inputFile:
        lines = inputFile.readlines()
        
        nrOfMolecules = int(lines[0].split()[1])
    
        molecules = []
        for i in range(1,nrOfMolecules+1):
            molecules.append(list(map(int,lines[i].split())))
        
        # create boolean matrix to indicate which atoms are part of same molecule:
        # TODO: Probably can still be done more efficiently?
        notInSameMolecule = np.ones((len(types), len(types)), dtype=bool)
        for i,mol in enumerate(molecules):
            for j,at in enumerate(mol):
                for k,at2 in enumerate(mol):
                    # print(at, at2)
                    notInSameMolecule[at,at2] = False
        
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
        
        nrOfDihedrals = int(lines[nrOfMolecules+nrOfBonds+nrOfAngles+3].split()[1])
        dihedrals = []
        dihedralConstants = []
        for i in range(nrOfMolecules+nrOfBonds+nrOfAngles+4, nrOfMolecules+nrOfBonds+nrOfAngles+nrOfDihedrals+4):
            dihedrals.append([int(lines[i].split()[0]),int(lines[i].split()[1]),int(lines[i].split()[2]),int(lines[i].split()[3])])
            dihedralConstants.append([float(lines[i].split()[4]),float(lines[i].split()[5]),float(lines[i].split()[6]),float(lines[i].split()[7])])
        dihedrals = np.asarray(dihedrals).reshape(nrOfDihedrals,4)
        dihedralConstants = np.asarray(dihedralConstants).reshape(nrOfDihedrals,4)
        
        sigma = []
        epsilon = []
        for i in range(nrOfMolecules+nrOfBonds+nrOfAngles+nrOfDihedrals+5, len(lines)):
            sigma.append(float(lines[i].split()[0]))
            epsilon.append(float(lines[i].split()[1]))
        sigma = np.asarray(sigma)
        epsilon = np.asarray(epsilon)
        
        return(molecules, notInSameMolecule, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon)

### DISTANCES & PBC ###

def distAtoms(positions):
    """ Computes distances between all atoms """
    diff = positions - positions[:,np.newaxis]
    dist = np.linalg.norm(diff,axis = 2)
    return(dist)

def distAtomsPBC(x):
    """ Computes distances between all atoms in closest copies, taking boundaries into account"""    
    diff = x - x[:,np.newaxis] 
    diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 

    # idea: if dist > 0.5*boxsize in some direction (x, y or z), then there is a closer copy. 
    # subtracting 0.5*boxsize in every direction where it is too large yields direction vector to closest neighbour
    # Still check correctness!
    dist = np.linalg.norm(diff,axis = 2)
    return(diff,dist) 
# Direction and distance are usually both needed, right? 
# Could also just return difference vector and do distance calculation elsewhere
#TODO replace difference vectors elsewhere, not only LJ?

def projectMolecules(x):
    """Projects entire molecules into box of given size
    
    Note: single atoms are not projected without the rest of their molecule
    Now computes weighted average for centre of mass. Is this worth the computation time?
    Might be that approximate centre of mass (e.g. unweighted) is good enough"""
    
    centers = np.zeros([len(types),3]) 
    for i in range(0, len(molecules)):
        centers[molecules[i]] = sum(x[molecules[i]] * m[molecules[i], np.newaxis] / sum(m[molecules[i]]))
    centersProj = centers % distAtomsPBC.boxSize
    x = x + (centersProj - centers)
    return x
    


### FORCES ###

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

# DIHEDRAL
def Vdihedral(t, C1, C2, C3, C4):
    """ Calculates periodic dihedral potential for given forcefield constants C """
    psi = np.pi - t
    return 0.5*(C1*(1+np.cos(psi)) + C2*(1-np.cos(2*psi)) + C3*(1+np.cos(3*psi)) + C4*(1-np.cos(4*psi)))

def Fdihedral(t, C1, C2, C3, C4):
    """ Calculates periodic dihedral force magnitude """
    psi = np.pi - t
    return 0.5*(-C1*np.sin(psi) + 2*C2*np.sin(2*psi) - 3*C3*np.sin(3*psi) + 4*C4*np.sin(4*psi))


### INTEGRATORS ###
def integratorEuler(x, v, a):
    """ Implementation of a single step for Euler integrator. """ 
    x = x + dt*v + (dt**2)/2*a
    v = v + dt*a
    return(x, v, a)

def integratorVerlocity(x, v, a):
    """ Implementation of a single step for Velocty Verlet integrator. """ 
    x_new = x + v*dt + (dt**2)/2*a
    a_new = computeForces(x_new, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon)/m[:,np.newaxis]
    v = v + dt/2*(a_new +a)
    return(x_new, v, a_new)

def integratorRK4(x, v, a):
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

# TODO: Fix integrators to work on only the input: x,v,a



### FORCES ###


def computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon):
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
        normalVec1 = np.cross(atomLeft-atomMiddle,np.cross(atomLeft-atomMiddle,atomRight-atomMiddle))
        normalVec2 = np.cross(atomMiddle-atomRight,np.cross(atomLeft-atomMiddle,atomRight-atomMiddle))
        
        # t = np.arccos(np.sum(dif1*dif2, axis = 1)/(np.linalg.norm(dif1, axis = 1)*np.linalg.norm(dif2, axis = 1))) # not corrected for floating point division
        t = np.arctan2(np.linalg.norm(np.cross(dif1,dif2), axis = 1),np.sum(dif1*dif2, axis = 1))
        
        Fangles = Fangle(t, angleConstants[:,0], angleConstants[:,1])
        
        FangleAtomLeft = Fangles[:,np.newaxis]/np.linalg.norm(atomLeft-atomMiddle, axis = 1)[:,np.newaxis] * normalVec1/np.linalg.norm(normalVec1, axis = 1)[:,np.newaxis] 
        FangleAtomRight = Fangles[:,np.newaxis]/np.linalg.norm(atomRight-atomMiddle, axis = 1)[:,np.newaxis] * normalVec2/np.linalg.norm(normalVec2, axis = 1)[:,np.newaxis]
        FangleAtomMiddle = - FangleAtomLeft - FangleAtomRight
        
        np.add.at(forces, angles[:,0], FangleAtomLeft)
        np.add.at(forces, angles[:,1], FangleAtomMiddle)
        np.add.at(forces, angles[:,2], FangleAtomRight)   
        
    # dihedrals
    if dihedrals.size > 0:
        rij = x[dihedrals[:,1]] - x[dihedrals[:,0]]
        rjk = x[dihedrals[:,2]] - x[dihedrals[:,1]] # TODO check right direction of these vectors and forces !!
        rkl = x[dihedrals[:,3]] - x[dihedrals[:,2]]
        ro = 0.5*x[dihedrals[:,1]] + 0.5*x[dihedrals[:,2]]
        rok = x[dihedrals[:,2]] - ro
        
        normalVec1 = np.cross(rij,rjk)
        normalVec2 = np.cross(rjk,rkl)
        
        # theta = np.arccos(np.sum(normalVec1*normalVec2, axis = 1)/(np.linalg.norm(normalVec1, axis = 1)*np.linalg.norm(normalVec2, axis = 1))) # not corrected for floating point division
        theta = np.arctan2(np.linalg.norm(np.cross(normalVec1,normalVec2), axis = 1),np.sum(normalVec1*normalVec2, axis = 1))*np.sign(np.sum(rij*normalVec2, axis = 1)) # corrected for floating point division
        
        Fdihedrals = Fdihedral(theta, dihedralConstants[:,0], dihedralConstants[:,1], dihedralConstants[:,2], dihedralConstants[:,3])
        angleAtomsijk = np.arccos(np.sum(rij*rjk, axis = 1)/(np.linalg.norm(rij, axis = 1)*np.linalg.norm(rjk, axis = 1))) # not corrected for floating point division
        angleAtomsjkl = np.arccos(np.sum(rjk*rkl, axis = 1)/(np.linalg.norm(rjk, axis = 1)*np.linalg.norm(rkl, axis = 1)))
        
        # angleAtomsijk = np.arctan2(np.linalg.norm(np.cross(rij,rjk), axis = 1),np.sum(rij*rjk, axis = 1))
        # angleAtomsjkl = np.arctan2(np.linalg.norm(np.cross(-rjk,rkl), axis = 1),np.sum(-rjk*rkl, axis = 1))
        
        FdihedralAtomi = -Fdihedrals[:,np.newaxis]/((np.linalg.norm(rij, axis = 1)*np.sin(angleAtomsijk))[:,np.newaxis]) * -normalVec1/np.linalg.norm(normalVec1, axis = 1)[:,np.newaxis]
        FdihedralAtoml = -Fdihedrals[:,np.newaxis]/((np.linalg.norm(rkl, axis = 1)*np.sin(angleAtomsjkl))[:,np.newaxis]) * normalVec2/np.linalg.norm(normalVec2, axis = 1)[:,np.newaxis] 
        
        FdihedralAtomk = (1/np.linalg.norm(rok, axis = 1)[:,np.newaxis]**2) * np.cross(-(np.cross(rok,FdihedralAtoml) + 0.5*np.cross(rkl, FdihedralAtoml) + 0.5*np.cross(rij, -FdihedralAtomi)),rok)
        
        FdihedralAtomj = - FdihedralAtomi - FdihedralAtomk - FdihedralAtoml
        # FdihedralAtomj = -FdihedralAtomi + np.sum(dif1*difCommon, axis = 1)/np.linalg.norm(difCommon, axis = 1)[:,np.newaxis]*FdihedralAtomi - np.sum(dif2*difCommon, axis = 1)/np.linalg.norm(dif2, axis = 1)[:,np.newaxis]*FdihedralAtoml
        # FdihedralAtomk = -FdihedralAtoml - np.sum(dif1*difCommon, axis = 1)/np.linalg.norm(difCommon, axis = 1)[:,np.newaxis]*FdihedralAtomi + np.sum(dif2*difCommon, axis = 1)/np.linalg.norm(dif2, axis = 1)[:,np.newaxis]*FdihedralAtoml
        
        np.add.at(forces, dihedrals[:,0], FdihedralAtomi)
        np.add.at(forces, dihedrals[:,1], FdihedralAtomj)
        np.add.at(forces, dihedrals[:,2], FdihedralAtomk)
        np.add.at(forces, dihedrals[:,3], FdihedralAtoml)
        
         # check that sum of torques is (approx) 0
        if checkForces:
            torquei = np.cross(x[dihedrals[:,0]] - ro, forces[dihedrals[:,0]])
            torquej = np.cross(x[dihedrals[:,1]] - ro, forces[dihedrals[:,1]])
            torquek = np.cross(x[dihedrals[:,2]] - ro, forces[dihedrals[:,2]])
            torquel = np.cross(x[dihedrals[:,3]] - ro, forces[dihedrals[:,3]])
            
            torqueSum = torquei + torquej + torquek + torquel
            
            if np.abs(np.sum(np.sum(torqueSum, axis = 0), axis = 0)) > 10**(-5): 
                print(f"Warning: sum of torques not equal to 0 but {np.sum(np.sum(torqueSum, axis = 0), axis = 0)}")
        
        
    # Lennard Jones forces
    if sigma.size > 0:
        # TODO update with PBC
        diff,dist = distAtomsPBC(x)
        # dist = distAtoms(x)
        # print(dist)
        
        # U = np.zeros((len(types), 3))
    
        e = np.sqrt(epsilon*epsilon[:,np.newaxis])
        s = 0.5*(sigma+sigma[:,np.newaxis])
        
        frac = np.divide(s, dist, out=np.zeros_like(s), where=dist!=0) # avoid division by 0
        frac6 = frac**6
        frac12 = frac6 ** 2
        U = 4*e*(frac12 - frac6) 
        V = np.sign(U)*np.divide(4*e*(6*frac6 - 12*frac12), dist, out=np.zeros_like(s), where=dist!=0) # avoid division by 0. 
        # TODO but what happens if dist=0? -> no division. That is ok since there is no forces on atoms at distance 0 anyway which is taken care of by notInSameMolecules
        # TODO is the sign correct? 
       
        #L = x-x[:,np.newaxis]
        L = diff
        
        # U = U*notInSameMolecule # these forces do not apply on atoms in the same molecule!
        # U = np.repeat(U, 3).reshape(len(types), len(types), 3)
        
        V = V*notInSameMolecule # these forces do not apply on atoms in the same molecule!
        V = np.repeat(V, 3).reshape(len(types), len(types), 3)
        
        #forces += np.sum(U*L, axis = 1)
        forces += np.sum(V*L, axis = 1)
    
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
    
    # forces check
    if checkForces:
        # check that sum of all forces is (approx) 0
        if np.abs(np.sum(np.sum(forces, axis = 0), axis = 0)) > 10**(-5): 
            print(f"Warning: sum of forces not equal to 0 but {np.sum(np.sum(forces, axis = 0), axis = 0)}")
    
    return(forces)

