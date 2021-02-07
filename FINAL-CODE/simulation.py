# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)

"""
import time as timer
import numpy as np
# import warnings
# warnings.filterwarnings("error") # treats warning as error; to stop computation when dividing by 0

### READ INPUT ###

# XYZ
def readXYZfile(fileName, timeStep): 
    """Read a .xyz file.
    
    Reads entire file for arbitrary number of timesteps.
    INPUT: .xyz file.
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
    massesDict = {'H': 1.0080, 'O': 15.9994, 'C': 12.0110} # in AMU (or Dalton )
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
        notInSameMolecule = np.ones((len(types), len(types)), dtype=bool)
        for i,mol in enumerate(molecules):
            for j,at in enumerate(mol):
                for k,at2 in enumerate(mol):
                    notInSameMolecule[at,at2] = False

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

def distAtomsPBC(x):
    """ Computes distances between all atoms in closest copies, taking boundaries into account."""    
    diff = x - x[:,np.newaxis] 
    diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 
    dist = np.linalg.norm(diff,axis = 2)
    return(diff,dist) 

def projectMolecules(x):
    """Projects entire molecules into box of given size
    
    Note: single atoms are not projected without the rest of their molecule."""
    global centers
    centers = np.zeros([len(types),3]) 
    for i in range(0, len(molecules)):
        centers[molecules[i]] = np.sum(x[molecules[i]], axis=0) / len(molecules[i]) # unweighted avg
    centersProj = centers % distAtomsPBC.boxSize
    x = x + (centersProj - centers)
    return x
    
### FORCES ###

# BOND
def Vbond(r, k, r0):
    """ Calculates harmonic bond potential."""
    return(1/2*k*(r-r0)**2)

def Fbond(r, k, r0):
    """ Calculates bond force magnitude."""
    return(-k*(r-r0))

# ANGLE
def Vangle(t, kt, t0):
    """ Calculates harmonic angular potential."""
    return 1/2*kt*(t-t0)**2

def Fangle(t, kt, t0):
    """ Calculates angular force magnitude."""
    return -kt*(t-t0)

# DIHEDRAL
def Vdihedral(t, C1, C2, C3, C4):
    """ Calculates periodic dihedral potential for given forcefield constants C."""
    psi = t - np.pi
    return 0.5*(C1*(1+np.cos(psi)) + C2*(1-np.cos(2*psi)) + C3*(1+np.cos(3*psi)) + C4*(1-np.cos(4*psi)))

def Fdihedral(t, C1, C2, C3, C4):
    """ Calculates periodic dihedral force magnitude."""
    psi =  t - np.pi
    return 0.5*(-C1*np.sin(psi) + 2*C2*np.sin(2*psi) - 3*C3*np.sin(3*psi) + 4*C4*np.sin(4*psi))

### INTEGRATORS ###
def integratorEuler(x, v):
    """ Implementation of a single step for Euler integrator.""" 
    forces, potential = computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigPair, epsPair, LJcutoff)
    a = forces/m[:,np.newaxis]
    x = x + dt*v + (dt**2)/2*a
    v = v + dt*a
    return(x, v, potential)

def integratorVerlet(x, xold, v):
    """ Implementation of a single step for Verlet integrator.""" 
    forces, potential = computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigPair, epsPair, LJcutoff)
    a = forces/m[:,np.newaxis]
    x_new = 2*x - xold + (dt**2)/2*a
    v = 1/(2*dt)*(x_new + xold)
    return(x_new, x, v, potential)

def integratorVerlocity(x, v, a):  
    """ Implementation of a single step for Velocty Verlet integrator.""" 
    x_new = x + v*dt + (dt**2)/2*a
    forces, potential = computeForces(x_new, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigPair, epsPair, LJcutoff)
    a_new = forces/m[:,np.newaxis]
    v = v + dt/2*(a_new + a)
    return(x_new, v, a_new, potential)

def integratorRK4(x, v):
    """ Implementation of a single step for Runge-Kutta order 4 integrator.""" 
    x1 = x + dt*v 
    f1, potential1 = computeForces(x1, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigPair, epsPair, LJcutoff)
    v1 = dt*f1/m[:,np.newaxis]
    x2 = x + dt/2*(v+v1/2)
    f2, potential2 = computeForces(x2, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigPair, epsPair, LJcutoff)
    v2 = dt*f2/m[:,np.newaxis]
    x3 = x + dt/2*(v+v2/2)
    f3, potential3 = computeForces(x3, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)
    v3 = dt*f3/m[:,np.newaxis]
    x4 = x + dt*(v+v3)
    f4, potential4 = computeForces(x4, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)
    v4 = dt*f4/m[:,np.newaxis]
    
    x = x + dt*v + dt**2*(v1+v2+v3)/6
    v = v + (v1+2*v2+2*v3+v4)/6
    
    f, potential = computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)
    
    return(x, v, potential1)

### FORCES ###

def computeBondForces(x, bonds, bondConstants):
    '''Computes forces and potentials acting on all bonds
    
    INPUT: Coordinates of all atoms, indices which atoms have bonds, bond force constants
    OUTPUT: bond forces acting on all particles, total bond potential over all particles
    '''
    r = np.linalg.norm(x[bonds[:,0]] - x[bonds[:,1]], axis = 1)
    Fbonds = Fbond(r, bondConstants[:,0], bondConstants[:,1])[:,np.newaxis]*(x[bonds[:,0]]-x[bonds[:,1]])/r[:,np.newaxis]
    
    pBonds = np.sum(Vbond(r,bondConstants[:,0], bondConstants[:,1]))
    return(Fbonds, pBonds)

def computeAngleForces(x, angles, angleConstants):
    '''Computes all angular forces and potentials
    
    INPUT: Coordinates of all atoms, indices which atoms form angles, angle force constants
    OUTPUT: angular forces acting on all particles, total angular potential over all particles
    '''
    atomLeft = x[angles[:,0]]
    atomMiddle = x[angles[:,1]]
    atomRight = x[angles[:,2]]
    
    dif1 = atomLeft-atomMiddle
    dif2 = atomRight-atomMiddle
    normalVec1 = np.cross(atomLeft-atomMiddle,np.cross(atomLeft-atomMiddle,atomRight-atomMiddle))
    normalVec2 = np.cross(atomMiddle-atomRight,np.cross(atomLeft-atomMiddle,atomRight-atomMiddle))
    
    t = np.arctan2(np.linalg.norm(np.cross(dif1,dif2), axis = 1),np.sum(dif1*dif2, axis = 1))
    
    Fangles = Fangle(t, angleConstants[:,0], angleConstants[:,1])
    
    FangleAtomLeft = Fangles[:,np.newaxis]/np.linalg.norm(atomLeft-atomMiddle, axis = 1)[:,np.newaxis] * normalVec1/np.linalg.norm(normalVec1, axis = 1)[:,np.newaxis] 
    FangleAtomRight = Fangles[:,np.newaxis]/np.linalg.norm(atomRight-atomMiddle, axis = 1)[:,np.newaxis] * normalVec2/np.linalg.norm(normalVec2, axis = 1)[:,np.newaxis]
    FangleAtomMiddle = - FangleAtomLeft - FangleAtomRight
    
    pAngles = np.sum(Vangle(t, angleConstants[:,0], angleConstants[:,1]))
    return(FangleAtomLeft, FangleAtomMiddle, FangleAtomRight, pAngles)

def computeDihedralForces():
    '''Computes dihedral forces and potentials
    
    INPUT: Coordinates of all atoms, indices which atoms form dihedrals, dihedral force constants
    OUTPUT: dihedral forces acting on all particles, total dihedral potential over all particles
    '''
    rij = x[dihedrals[:,1]] - x[dihedrals[:,0]]
    rjk = x[dihedrals[:,2]] - x[dihedrals[:,1]] 
    rkl = x[dihedrals[:,3]] - x[dihedrals[:,2]]
    ro = 0.5*x[dihedrals[:,1]] + 0.5*x[dihedrals[:,2]]
    rok = x[dihedrals[:,2]] - ro
    
    normalVec1 = np.cross(rij,rjk)
    normalVec2 = np.cross(rjk,rkl)
    
    theta = np.arctan2(np.linalg.norm(np.cross(normalVec1,normalVec2), axis = 1),np.sum(normalVec1*normalVec2, axis = 1))*np.sign(np.sum(rij*normalVec2, axis = 1)) # corrected for floating point division
    
    Fdihedrals = Fdihedral(theta, dihedralConstants[:,0], dihedralConstants[:,1], dihedralConstants[:,2], dihedralConstants[:,3])
    angleAtomsijk = np.arctan2(np.linalg.norm(np.cross(rij,rjk), axis = 1),np.sum(rij*rjk, axis = 1)) #wrong sign?
    angleAtomsjkl = np.arctan2(np.linalg.norm(np.cross(rjk,rkl), axis = 1),np.sum(rjk*rkl, axis = 1))
    
    FdihedralAtomi = Fdihedrals[:,np.newaxis]/((np.linalg.norm(rij, axis = 1)*np.sin(angleAtomsijk))[:,np.newaxis]) * -normalVec1/np.linalg.norm(normalVec1, axis = 1)[:,np.newaxis]
    FdihedralAtoml = Fdihedrals[:,np.newaxis]/((np.linalg.norm(rkl, axis = 1)*np.sin(angleAtomsjkl))[:,np.newaxis]) * normalVec2/np.linalg.norm(normalVec2, axis = 1)[:,np.newaxis] 
    
    FdihedralAtomk = (1/np.linalg.norm(rok, axis = 1)[:,np.newaxis]**2) * np.cross(-(np.cross(rok,FdihedralAtoml) + 0.5*np.cross(rkl, FdihedralAtoml) + 0.5*np.cross(rij, -FdihedralAtomi)),rok)
    
    FdihedralAtomj = - FdihedralAtomi - FdihedralAtomk - FdihedralAtoml
    
    pDihedrals = -np.sum(Vdihedral(theta, dihedralConstants[:,0], dihedralConstants[:,1], dihedralConstants[:,2], dihedralConstants[:,3])) 
    return(FdihedralAtomi, FdihedralAtomj, FdihedralAtomk, FdihedralAtoml, pDihedrals)

def computeLJForces(x, sigma, epsilon, LJcutoff):
    '''Computes Lennard-Jones forces and potentials for all pairs closer together than given cutoff
    
    INPUT: Coordinates of all atoms, pairwise mixed LJ parameter arrays (precomputed), cutoff length
    OUTPUT: LJ forces acting on all particles, total LJ potential over all particles, both considering cutoff
    '''
    diff,dist = distAtomsPBC(x)
    atomPairs = np.where(np.multiply(dist < LJcutoff, np.triu(notInSameMolecule)) == True) #tuple of arrays; sort of adjacency list. triu to avoid duplicates
    distReciprocal = np.power(dist[atomPairs[0],atomPairs[1]], -1)
    frac = np.multiply(sigPair[atomPairs[0],atomPairs[1]], distReciprocal) #sigPair is precomputed for all pairs of sigma
    frac6 = np.power(frac, 6)
    e = epsPair[atomPairs[0],atomPairs[1]] #epsPair is precomputed for all pairs of epsilons
    
    epsFrac6 = np.multiply(4*e,frac6) # to re-use in computation of both U and V
    epsFrac12 = np.multiply(epsFrac6, frac6)
    
    U = epsFrac12 - epsFrac6 # = 4*e*(frac12 - frac6); potential
    
    V = np.multiply((12*epsFrac12 - 6*epsFrac6), distReciprocal) # = - (4*e*(6*frac6 - 12*frac12) / dist
    
    LJforces = np.multiply(np.multiply(diff[atomPairs[0],atomPairs[1]], distReciprocal[:,np.newaxis]), V[:,np.newaxis])
    
    pLJ = np.sum(U) 
    return(LJforces, pLJ, atomPairs)

def computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff):
    """Calculates all forces acting on all particles using force functions and parameters from topology file."""
    forces = np.zeros((len(types),3), dtype = float)
    potentials = np.zeros(4, dtype = float)
    
    # bonds
    if bonds.size > 0:
        Fbonds, pBonds = computeBondForces(x, bonds, bondConstants)
        np.add.at(forces, bonds[:,0], Fbonds)
        np.add.at(forces, bonds[:,1], -Fbonds)  
        potentials[0] = pBonds
    
    # angles 
    if angles.size > 0:
        FangleAtomLeft, FangleAtomMiddle, FangleAtomRight, pAngles = computeAngleForces(x, angles, angleConstants)
        np.add.at(forces, angles[:,0], FangleAtomLeft)
        np.add.at(forces, angles[:,1], FangleAtomMiddle)
        np.add.at(forces, angles[:,2], FangleAtomRight)  
        potentials[1] = pAngles
        
    # dihedrals
    if dihedrals.size > 0:
        FdihedralAtomi, FdihedralAtomj, FdihedralAtomk, FdihedralAtoml, pDihedrals = computeDihedralForces()
        np.add.at(forces, dihedrals[:,0], FdihedralAtomi)
        np.add.at(forces, dihedrals[:,1], FdihedralAtomj)
        np.add.at(forces, dihedrals[:,2], FdihedralAtomk)
        np.add.at(forces, dihedrals[:,3], FdihedralAtoml)
        potentials[2] = pDihedrals
        
    # Lennard Jones forces
    if sigPair.size > 0:  
        LJforces, pLJ, atomPairs = computeLJForces(x, sigma, epsilon, LJcutoff)
        np.add.at(forces, atomPairs[0], -LJforces) 
        np.add.at(forces, atomPairs[1], LJforces)
        potentials[3] = pLJ
    
    return(forces, potentials)


### PARAMETERS ###
sizesSmall = {'Water': 31.08, 'Ethanol': 32.22, 'Mixture': 32.29} #Angstrom
sizesLarge = {'Water': 49.72, 'Ethanol': 50.61, 'Mixture': 51.65} #Angstrom
# boxsizes are chosen so that all grid positions are filled
                
def setSimulation(substance, small = True, therm = True):
    '''Sets parameters for different simulations
    
    INPUT: substance (water/ethanol/mixture), box size (small/large), whether to use a thermostat (false/true)
    OUTPUT: None; global variables set: boxsize for PBC, names of input files, name for output file, boolean 'thermostat'
    '''
    global inputTimeStep, inputFileName, topologyFileName, outputFileName, thermostat
    
    if small:
        distAtomsPBC.boxSize = sizesSmall.get(substance)
    else: 
        distAtomsPBC.boxSize = sizesLarge.get(substance)
    
    if therm:
        thermo = 'Thermostat'
    else:
        thermo = ''
    
    size = str(distAtomsPBC.boxSize)
    inputTimeStep = 0
    inputFileName = substance + size + 'Initial.xyz'
    topologyFileName = substance + size + 'Topology.txt'
    outputFileName = substance + size + thermo 
    thermostat = therm

setSimulation('Ethanol') # choose 'Water", 'Ethanol' or 'Mixture'

### SIMULATION ###
types, x, m = readXYZfile(inputFileName, inputTimeStep)
molecules, notInSameMolecule, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon = readTopologyFile(topologyFileName)
LJcutoff = 2.5*np.max(sigma) #Angstrom; advised in literature: 2.5*sigma

time = 0 # ps
endTime = 1 # ps; should be 1ns = 1000ps in final simulation or 0.1ns = 100ps 
dt = 0.002 # ps; suggestion was to start at 2fs for final simulations, paper uses 0.5fs

u = np.random.uniform(size=3*len(types)).reshape((len(types),3)) # random starting velocity vector
u = u/np.linalg.norm(u,axis = 1)[:,np.newaxis] # normalize
v = 0.01*u # A/ps

# For Gaussian Thermostat:
temperatureDesired = 298.15 # Kelvin 
kB = 0.8314459727525677 #[A^2 AMU] / [ps^2 K]; from regular kB via 1.38064852  * 6.02214076 * 10 **(-23 + 20 - 24 + 26) 
Nf = 3*len(x) # 3*, as atoms have 3D velocity vector and only translational freedom matters
c = 1/( kB * Nf) 
    
# For measuring:
nrSteps = int(np.ceil(endTime/dt))
Ekin = np.zeros(nrSteps)
Epot = np.zeros((nrSteps,4))
temperatures = np.zeros(nrSteps)
j = 0

# Precomputations: 
sigPair = 0.5*(sigma + sigma[:, np.newaxis]) 
epsPair = np.sqrt(epsilon * epsilon[:, np.newaxis])
f_init, pot = computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigPair, epsPair, LJcutoff)
a = f_init / m[:,np.newaxis]
temperatureSystem = np.sum(m * np.linalg.norm(v, axis = 1)**2) / (Nf * kB)
v = v * np.sqrt(temperatureDesired/temperatureSystem) 

#Simulation loop:
with open(outputFileName + 'Output.xyz', "w") as outputFile: # clear file
    outputFile.write("") 
with open(outputFileName + 'Measurables.txt', "w") as outputFileMeas: # clear file
    outputFileMeas.write("") 
simStartTime = timer.time()
with open(outputFileName + 'Output.xyz', "a") as outputFile:
    with open(outputFileName + 'Measurables.txt', "a") as measurables:
        measurables.write("Kinetic,  Bond potential, Angle pot, Dihedral pot, LJ pot, Temp. before scaling, temperature after \n")
    
        while (time < endTime) : 
            print(time, " out of ", endTime)
            if (time % (30*dt) < dt): #to print every 30th frame. '==0' does not work, as floats are not exact. add 'or True' to print all
                outputFile.write(f"{len(types)}\n")
                outputFile.write(f"This is a comment and the time is {time:5.4f}\n")
                for i, atom in enumerate(x):
                    outputFile.write(f"{types[i]} {x[i,0]:10.5f} {x[i,1]:10.5f} {x[i,2]:10.5f}\n")  
            
            
            x, v, a, potentials = integratorVerlocity(x, v, a) #integration
            x = projectMolecules(x) #projection
            time += dt
            
            temperatureSystem = np.sum(m * np.linalg.norm(v, axis = 1)**2) / (Nf * kB)
            if thermostat: #velocity rescaling
                v = v * np.sqrt(temperatureDesired/temperatureSystem) 
            
            # measurables
            EkinSyst = np.sum(0.5 * m * (np.linalg.norm(v, axis=1)**2))
            tempAfter = 2*EkinSyst / (Nf * kB)
            Ekin[j] = EkinSyst
            Epot[j] = potentials
            temperatures[j] = temperatureSystem
            measurables.write(f"{EkinSyst:10.5f} {potentials[0]:10.5f} {potentials[1]:10.5f} {potentials[2]:10.5f} {potentials[3]:10.5f} {temperatureSystem:10.5f} {tempAfter:10.5f} \n")
            j += 1
        
duration = timer.time() - simStartTime
print("Simulation duration was ", int(duration/3600), 'hours, ', int((duration%3600)/60), " minutes and ", int(duration%60), "seconds")  