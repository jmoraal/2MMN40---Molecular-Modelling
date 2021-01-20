# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 21:17:16 2021

@author: s161981
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@authors: Sanne van Kempen (1017389) & Jan Moraal (1016866)

"""
import time as timer
import numpy as np
# import warnings
# warnings.filterwarnings("error") # allows try/except-constructions for warnings (i.e. stop computation when dividing by 0)
#import functions as fun # idea was to put all fcts there, but too many non-parametrised variables


### READ INPUT ###

# XYZ:
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
    """ Computes distances between all atoms in closest copies, taking boundaries into account"""    
    diff = x - x[:,np.newaxis] 
    diff = diff - np.floor(0.5 + diff/distAtomsPBC.boxSize)*distAtomsPBC.boxSize 

    dist = np.linalg.norm(diff,axis = 2)
    return(diff,dist) 


def projectMolecules(x):
    """Projects entire molecules into box of given size
    
    Note: single atoms are not projected without the rest of their molecule
    Now computes weighted average for centre of mass. Is this worth the computation time?
    Might be that approximate centre of mass (e.g. unweighted) is good enough"""
    #TODO make faster?
    centers = np.zeros([len(types),3]) 
    #centers[molecules] = sum(x[molecules]  / sum([molecules])) #should be doable without forloop
    for i in range(0, len(molecules)):
        # centers[molecules[i]] = sum(x[molecules[i]] * m[molecules[i], np.newaxis]) / sum(m[molecules[i]]) # weighted avg
        centers[molecules[i]] = np.sum(x[molecules[i]], axis=0) / len(molecules[i]) # unweighted avg
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
    return(x, v) #a was previously returned too, but why?

# def integratorVerlocity(x, v, a, a_new):  #TODO: should be able to do this with these parameters w/o computeForces for cleaner code
#     """ Implementation of a single step for Velocty Verlet integrator. """ 
#     x_new = x + v*dt + (dt**2)/2*a
#     #a_new = computeForces(x_new, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)/m[:,np.newaxis]
#     v = v + dt/2*(a_new + a)
#     return(x_new, v)

def integratorVerlocity(x, v, a):  
    """ Implementation of a single step for Velocty Verlet integrator. """ 
    x_new = x + v*dt + (dt**2)/2*a
    forces, potential = computeForces(x_new, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)
    a_new = forces/m[:,np.newaxis]
    v = v + dt/2*(a_new + a)
    return(x_new, v, a_new, potential)

def integratorRK4(x, v, a):
    """ Implementation of a single step for Runge-Kutta order 4 integrator. """ 
    x1 = x + dt*v + (dt**2)/2*a 
    v1 = dt*computeForces(x1, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)/m[:,np.newaxis]
    x2 = x + dt/2*(v+v1/2) + (dt**2)/2*a 
    v2 = dt*computeForces(x2, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)/m[:,np.newaxis]
    x3 = x + dt/2*(v+v2/2) + (dt**2)/2*a
    v3 = dt*computeForces(x3, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)/m[:,np.newaxis]
    x4 = x + dt*(v+v3) + (dt**2)/2*a 
    v4 = dt*computeForces(x4, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)/m[:,np.newaxis]
    
    v = v + (v1+2*v2+2*v3+v4)/6
    return(x1, v, a)



### FORCES ###


def computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff):
    """Caltulate forces in one go with help of topology file."""
    forces = np.zeros((len(types),3), dtype = float)
    potentials = np.zeros(4, dtype = float)
    
    # bonds
    if bonds.size > 0:
        r = np.linalg.norm(x[bonds[:,0]] - x[bonds[:,1]], axis = 1)
        Fbonds = Fbond(r, bondConstants[:,0], bondConstants[:,1])[:,np.newaxis]*(x[bonds[:,0]]-x[bonds[:,1]])/r[:,np.newaxis]
        np.add.at(forces, bonds[:,0], Fbonds)
        np.add.at(forces, bonds[:,1], -Fbonds)  
        
        potentials[0] = np.sum(Vbond(r,bondConstants[:,0], bondConstants[:,1]))
    
    # angles 
    if angles.size > 0:
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
        
        np.add.at(forces, angles[:,0], FangleAtomLeft)
        np.add.at(forces, angles[:,1], FangleAtomMiddle)
        np.add.at(forces, angles[:,2], FangleAtomRight)   
        
        potentials[1] = np.sum(Vangle(t, angleConstants[:,0], angleConstants[:,1]))
        
    # dihedrals
    if dihedrals.size > 0:
        rij = x[dihedrals[:,1]] - x[dihedrals[:,0]]
        rjk = x[dihedrals[:,2]] - x[dihedrals[:,1]] 
        rkl = x[dihedrals[:,3]] - x[dihedrals[:,2]]
        ro = 0.5*x[dihedrals[:,1]] + 0.5*x[dihedrals[:,2]]
        rok = x[dihedrals[:,2]] - ro
        
        normalVec1 = np.cross(rij,rjk)
        normalVec2 = np.cross(rjk,rkl)
        
        theta = np.arctan2(np.linalg.norm(np.cross(normalVec1,normalVec2), axis = 1),np.sum(normalVec1*normalVec2, axis = 1))*np.sign(np.sum(rij*normalVec2, axis = 1)) # corrected for floating point division
        
        Fdihedrals = Fdihedral(theta, dihedralConstants[:,0], dihedralConstants[:,1], dihedralConstants[:,2], dihedralConstants[:,3])
        angleAtomsijk = np.arccos(np.sum(rij*rjk, axis = 1)/(np.linalg.norm(rij, axis = 1)*np.linalg.norm(rjk, axis = 1))) # not corrected for floating point division
        angleAtomsjkl = np.arccos(np.sum(rjk*rkl, axis = 1)/(np.linalg.norm(rjk, axis = 1)*np.linalg.norm(rkl, axis = 1)))
        
        FdihedralAtomi = -Fdihedrals[:,np.newaxis]/((np.linalg.norm(rij, axis = 1)*np.sin(angleAtomsijk))[:,np.newaxis]) * -normalVec1/np.linalg.norm(normalVec1, axis = 1)[:,np.newaxis]
        FdihedralAtoml = -Fdihedrals[:,np.newaxis]/((np.linalg.norm(rkl, axis = 1)*np.sin(angleAtomsjkl))[:,np.newaxis]) * normalVec2/np.linalg.norm(normalVec2, axis = 1)[:,np.newaxis] 
        
        FdihedralAtomk = (1/np.linalg.norm(rok, axis = 1)[:,np.newaxis]**2) * np.cross(-(np.cross(rok,FdihedralAtoml) + 0.5*np.cross(rkl, FdihedralAtoml) + 0.5*np.cross(rij, -FdihedralAtomi)),rok)
        
        FdihedralAtomj = - FdihedralAtomi - FdihedralAtomk - FdihedralAtoml
        
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
                
        potentials[2] = np.sum(Vdihedral(theta, dihedralConstants[:,0], dihedralConstants[:,1], dihedralConstants[:,2], dihedralConstants[:,3]))
        
        
    # Lennard Jones forces
    if sigma.size > 0:
        diff,dist = distAtomsPBC(x)
    
        e = np.sqrt(epsilon*epsilon[:,np.newaxis])
        s = 0.5*(sigma+sigma[:,np.newaxis])
        
        frac = np.divide(s, dist, out=np.zeros_like(s), where=dist!=0) # avoid division by 0
        frac6 = frac**6
        frac12 = frac6 ** 2
        U = 4*e*(frac12 - frac6) # potential
        V = np.sign(U)*np.divide(4*e*(6*frac6 - 12*frac12), dist, out=np.zeros_like(s), where=dist!=0) # avoid division by 0. 
        # TODO is the sign correct? 
       
        L = diff
        V = V*notInSameMolecule*(dist < LJcutoff) # these forces do not apply on atoms in the same molecule, and only apply when dist < cutoff
        V = np.repeat(V, 3).reshape(len(types), len(types), 3)
        
        forces += np.sum(V*L, axis = 1)
        
        potentials[3] = np.sum(U)
    
    
    # forces check
    if checkForces:
        # check that sum of all forces is (approx) 0
        if np.abs(np.sum(np.sum(forces, axis = 0), axis = 0)) > 10**(-5): 
            print(f"Warning: sum of forces not equal to 0 but {np.sum(np.sum(forces, axis = 0), axis = 0)}")
    
    return(forces, potentials)

sizesSmall = {'Water': 31.08, 'Ethanol': 32.22, 'Mixture': 32.29}
sizesLarge = {'Water': 49.72, 'Ethanol': 50.61, 'Mixture': 51.65}
# boxsizes are chosen so that all grid positions are filled
                
### PARAMETERS ###
def setSimulation(substance, small = True, therm = True):
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
    outputFileName = substance + size + thermo + 'Output.xyz'
    thermostat = therm

setSimulation('Water')

# inputFileName = "MixedMolecules.xyz"
# inputTimeStep = 0
# topologyFileName = "MixedMoleculesTopology.txt"
# outputFileName = "MixedMoleculesOutput.xyz"
# distAtomsPBC.boxSize = 48.42
# thermostat = False


### SIMULATION ###
types, x, m = readXYZfile(inputFileName, inputTimeStep)
molecules, notInSameMolecule, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon = readTopologyFile(topologyFileName)
LJcutoff = 2.5*np.max(sigma)

time = 0 #ps
endTime = 0.03 #ps; should be 1ns = 1000ps in final simulation
dt = 0.003 #ps; suggestion was to start at 2fs for final simulations, larger might be better (without exploding at least)

u = np.random.uniform(size=3*len(types)).reshape((len(types),3)) # random starting velocity vector
u = u/np.linalg.norm(u,axis = 1)[:,np.newaxis] # normalize
v = 0.01*u # A/ps

#For Gaussian Thermostat:
if thermostat: 
    temperatureDesired = 298.15 # Kelvin
    kB = 1.38064852 * 1.6605390666 # [A^2 AMU]/[K ps^2]; hence -20+23-27+24 = 0 'in the exponent'
    Nf = 6*len(x) # 6*, as atoms have 3D position and velocity vector so 6 degrees of freedom
    c = 2/(3*kB*len(x)) #for alternative computation using equipartition theorem
    
    
# For measuring:
Ekin = []
Epot = []
checkForces = False

with open(outputFileName, "w") as outputFile: # clear file
    outputFile.write("") 
simStartTime = timer.time()
with open(outputFileName, "a") as outputFile:
    f_init, pot = computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)
    a = f_init / m[:,np.newaxis]
    
    while (time <= endTime) : 
        print(time, " out of ", endTime)
        outputFile.write(f"{len(types)}\n")
        outputFile.write(f"This is a comment and the time is {time:5.4f}\n")
        for i, atom in enumerate(x):
            outputFile.write(f"{types[i]} {x[i,0]:10.5f} {x[i,1]:10.5f} {x[i,2]:10.5f}\n")  
        
        #forces, potentials = computeForces(x, bonds, bondConstants, angles, angleConstants, dihedrals, dihedralConstants, sigma, epsilon, LJcutoff)
        #a_new = forces / m[:,np.newaxis]
        #x, v = integratorVerlocity(x, v, a, a_new)
        #a = a_new
        
        x, v, a, potentials = integratorVerlocity(x, v, a)
        x = projectMolecules(x) #TODO is this the right place, or should it be before integration/force computation?
        time += dt
        
        if thermostat: 
            temperatureSystem = c*np.sum((np.linalg.norm(v, axis=1)**2)/m) #via equipartition theorem
            v = v * np.sqrt(temperatureDesired/temperatureSystem) 
        
        # measurables
        EkinSyst = np.sum(0.5 * m * (np.linalg.norm(v, axis=1)**2))
        Ekin.append(EkinSyst)
        EpotSyst =  np.sum(potentials)
        Epot.append(EpotSyst)
        

duration = timer.time() - simStartTime
print("Simulation duration was ", duration, " seconds")
         

### PLOT ENERGY ###

# import matplotlib.pyplot as plt
# t = np.arange(0,time-dt, dt)
# Etot = np.array(Ekin) + np.array(Epot)
# plt.plot(t, Ekin, label = 'Kinetic energy')
# plt.plot(t, Epot, label = 'Potential energy')
# plt.plot(t, Etot, label = 'Total energy')
# plt.title('Energy in the system')
# plt.xlabel('Time (ps)')
# plt.ylabel('Energy')
# plt.legend()
# plt.show()