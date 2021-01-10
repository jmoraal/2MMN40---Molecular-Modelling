# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:54:00 2021

@author: s161981
"""
import numpy as np

### Water ###
def writeWaterTopology(nrOfMolecules, boxSize): 
    outputFileName = 'Water' + str(boxSize) + 'Topology.txt'
    
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"molecules {nrOfMolecules}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f"{3*i} {3*i+1} {3*i+2}\n")
        
        outputFile.write(f"bonds {nrOfMolecules*2}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f"{3*i} {3*i+1} 5024.16 0.9572 \n")
            outputFile.write(f"{3*i} {3*i+2} 5024.16 0.9572 \n")
        
        outputFile.write(f"angles {nrOfMolecules}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f"{3*i+1} {3*i} {3*i+2} 628.02 1.8242181 \n")
            
        outputFile.write("dihedrals 0\n")
        
        outputFile.write(f"LJ {nrOfMolecules*3}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write("3.15061 0.66386\n")
            outputFile.write("0 0\n")
            outputFile.write("0 0\n")
  
    
### Ethanol ###
def writeEthanolTopology(nrOfMolecules, boxSize): 
    outputFileName = 'Ethanol' + str(boxSize) + 'Topology.txt'
    
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
        
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"molecules {nrOfMolecules}\n")
        for i in range(0,nrOfMolecules):
            for j in range(0,9):
                outputFile.write(f"{9*i + j} ")
            outputFile.write("\n")
        
        outputFile.write(f"bonds {nrOfMolecules*8}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f"{9*i +1} {9*i + 0} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +2} {9*i + 0} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +3} {9*i + 0} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +0} {9*i + 4} 2242.624 1.529 \n" )
            outputFile.write(f"{9*i +5} {9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +6} {9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +4} {9*i + 7} 2677.76 1.41 \n" )
            outputFile.write(f"{9*i +7} {9*i + 8} 4627.5 0.945 \n" )
        
        outputFile.write(f"angles {nrOfMolecules*13}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f'{13*i +1} {13*i + 0} {13*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{13*i +2} {13*i + 0} {13*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{13*i +3} {13*i + 0} {13*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{13*i +3} {13*i + 0} {13*i + 2} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +3} {13*i + 0} {13*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +2} {13*i + 0} {13*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +5} {13*i + 4} {13*i + 6} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +0} {13*i + 4} {13*i + 6} 313.8 1.93208 \n' )
            outputFile.write(f'{13*i +0} {13*i + 4} {13*i + 5} 313.8 1.93208 \n' )
            outputFile.write(f'{13*i +0} {13*i + 4} {13*i + 7} 414.4 1.91114 \n' )
            outputFile.write(f'{13*i +4} {13*i + 7} {13*i + 8} 460.24 1.89368 \n' )
            outputFile.write(f'{13*i +5} {13*i + 4} {13*i + 7} 292.88 1.91114 \n' )
            outputFile.write(f'{13*i +6} {13*i + 4} {13*i + 7} 292.88 1.91114 \n' )
    
        outputFile.write(f"dihedrals {nrOfMolecules*12}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f'{12*i +1} {12*i + 0} {12*i + 4} {12*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +2} {12*i + 0} {12*i + 4} {12*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +3} {12*i + 0} {12*i + 4} {12*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +1} {12*i + 0} {12*i + 4} {12*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +2} {12*i + 0} {12*i + 4} {12*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +3} {12*i + 0} {12*i + 4} {12*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +1} {12*i + 0} {12*i + 4} {12*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{12*i +2} {12*i + 0} {12*i + 4} {12*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{12*i +3} {12*i + 0} {12*i + 4} {12*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{12*i +0} {12*i + 4} {12*i + 7} {12*i + 8} -0.4431 3.83255 0.72801 -4.11705 \n' )
            outputFile.write(f'{12*i +5} {12*i + 4} {12*i + 7} {12*i + 8} 0.9414 2.8242 0 -3.7656 \n' )
            outputFile.write(f'{12*i +6} {12*i + 4} {12*i + 7} {12*i + 8} 0.9414 2.8242 0 -3.7656 \n' )
        
        outputFile.write(f"LJ {nrOfMolecules*9}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write('3.5 0.276144 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('3.5 0.276144 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('3.12 0.71128 \n')
            outputFile.write('0 0 \n')


### Mixture ###
def writeMixtureTopology(nrOfMolecules, boxSize): 
    outputFileName = 'Mixture' + str(boxSize) + 'Topology.txt'
    nrEthanol = int(nrOfMolecules * 228 / 3716) #based on 14.3% mass ethanol
    nrWater = nrOfMolecules - nrEthanol
    
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"molecules {nrWater + nrEthanol}\n")
        for i in range(0,nrWater):
            outputFile.write(f"{3*i} {3*i+1} {3*i+2}\n")
        for i in range(0,nrEthanol):
            for j in range(0,9):
                outputFile.write(f"{9*i + j} ")
            outputFile.write("\n")
        
        outputFile.write(f"bonds {nrWater*2 + nrEthanol * 8}\n")
        for i in range(0,nrWater):
            outputFile.write(f"{3*i} {3*i+1} 5024.16 0.9572 \n")
            outputFile.write(f"{3*i} {3*i+2} 5024.16 0.9572 \n")
        for i in range(0,nrEthanol):
            outputFile.write(f"{9*i +1} {9*i + 0} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +2} {9*i + 0} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +3} {9*i + 0} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +0} {9*i + 4} 2242.624 1.529 \n" )
            outputFile.write(f"{9*i +5} {9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +6} {9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +4} {9*i + 7} 2677.76 1.41 \n" )
            outputFile.write(f"{9*i +7} {9*i + 8} 4627.5 0.945 \n" )
        
        outputFile.write(f"angles {nrWater + nrEthanol * 13}\n")
        for i in range(0,nrWater):
            outputFile.write(f"{3*i+1} {3*i} {3*i+2} 628.02 1.8242181 \n")
        for i in range(0,nrEthanol):
            outputFile.write(f'{13*i +1} {13*i + 0} {13*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{13*i +2} {13*i + 0} {13*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{13*i +3} {13*i + 0} {13*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{13*i +3} {13*i + 0} {13*i + 2} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +3} {13*i + 0} {13*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +2} {13*i + 0} {13*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +5} {13*i + 4} {13*i + 6} 276.144 1.88146 \n' )
            outputFile.write(f'{13*i +0} {13*i + 4} {13*i + 6} 313.8 1.93208 \n' )
            outputFile.write(f'{13*i +0} {13*i + 4} {13*i + 5} 313.8 1.93208 \n' )
            outputFile.write(f'{13*i +0} {13*i + 4} {13*i + 7} 414.4 1.91114 \n' )
            outputFile.write(f'{13*i +4} {13*i + 7} {13*i + 8} 460.24 1.89368 \n' )
            outputFile.write(f'{13*i +5} {13*i + 4} {13*i + 7} 292.88 1.91114 \n' )
            outputFile.write(f'{13*i +6} {13*i + 4} {13*i + 7} 292.88 1.91114 \n' )
            
        outputFile.write(f"dihedrals {nrEthanol * 12}\n")
        for i in range(0,nrEthanol):
            outputFile.write(f'{12*i +1} {12*i + 0} {12*i + 4} {12*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +2} {12*i + 0} {12*i + 4} {12*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +3} {12*i + 0} {12*i + 4} {12*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +1} {12*i + 0} {12*i + 4} {12*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +2} {12*i + 0} {12*i + 4} {12*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +3} {12*i + 0} {12*i + 4} {12*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{12*i +1} {12*i + 0} {12*i + 4} {12*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{12*i +2} {12*i + 0} {12*i + 4} {12*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{12*i +3} {12*i + 0} {12*i + 4} {12*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{12*i +0} {12*i + 4} {12*i + 7} {12*i + 8} -0.4431 3.83255 0.72801 -4.11705 \n' )
            outputFile.write(f'{12*i +5} {12*i + 4} {12*i + 7} {12*i + 8} 0.9414 2.8242 0 -3.7656 \n' )
            outputFile.write(f'{12*i +6} {12*i + 4} {12*i + 7} {12*i + 8} 0.9414 2.8242 0 -3.7656 \n' )
        
        outputFile.write(f"LJ {nrWater*3 + nrEthanol * 9}\n")
        for i in range(0,nrWater):
            outputFile.write("3.15061 0.66386\n")
            outputFile.write("0 0\n")
            outputFile.write("0 0\n")
        for i in range(0,nrEthanol):
            outputFile.write('3.5 0.276144 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('3.5 0.276144 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('2.5 0.12552 \n')
            outputFile.write('3.12 0.71128 \n')
            outputFile.write('0 0 \n') 
            


### Generate xyz's ###
def writeWaterXYZ(nrOfMolecules, boxSize):
    outputFileName = 'Water' + str(boxSize) + 'Initial.xyz'
    d = boxSize / (nrOfMolecules)**(1/3)
    n = int(nrOfMolecules**(1/3)) + 2
    rand = []
    for i in range(1,n): 
        for j in range(1,n): 
            for k in range(1,n): 
                rand.append([(i+1)*d, (j+1)*d, (k+1)*d])
        
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"{nrOfMolecules*3}\n")
        outputFile.write("Comment t = 0\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f"O {rand[i][0]} {rand[i][1]} {rand[i][2]}\n")
            outputFile.write(f"H {rand[i][0]-1} {rand[i][1]} {rand[i][2]}\n")
            outputFile.write(f"H {rand[i][0]} {rand[i][1]-0.3} {rand[i][2]-0.9}\n")
            
def writeEthanolXYZ(nrOfMolecules, boxSize):
    outputFileName = 'Ethanol' + str(boxSize) + 'Initial.xyz'
    d = boxSize / (nrOfMolecules)**(1/3)
    n = int(nrOfMolecules**(1/3)) + 2
    rand = []
    for i in range(1,n): 
        for j in range(1,n): 
            for k in range(1,n): 
                rand.append([(i+1)*d, (j+1)*d, (k+1)*d])
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"{nrOfMolecules*9}\n")
        outputFile.write("Comment t = 0\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f'C {rand[i][0] + -1} {rand[i][1] + -1} {rand[i][2] + 0.5} \n')
            outputFile.write(f'H {rand[i][0] + -1} {rand[i][1] + -2} {rand[i][2] + 0.2} \n')
            outputFile.write(f'H {rand[i][0] + -2} {rand[i][1] + -1} {rand[i][2] + 1.1} \n')
            outputFile.write(f'H {rand[i][0] + -1} {rand[i][1] + 0} {rand[i][2] + 0.1} \n')
            outputFile.write(f'C {rand[i][0] + 0.5} {rand[i][1] + -1} {rand[i][2] + 0.9} \n')
            outputFile.write(f'H {rand[i][0] + 0.5} {rand[i][1] + -0.2} {rand[i][2] + 1.7} \n')
            outputFile.write(f'H {rand[i][0] + 0.5} {rand[i][1] + -1.7} {rand[i][2] + 1.7} \n')
            outputFile.write(f'O {rand[i][0] + 2} {rand[i][1] + -1} {rand[i][2] + 0.4} \n')
            outputFile.write(f'H {rand[i][0] + 2.5} {rand[i][1] + -0.5} {rand[i][2] + 0.2} \n')


def writeMixtureXYZ(nrOfMolecules, boxSize):
    outputFileName = 'Mixture' + str(boxSize) + 'Initial.xyz'
    nrEthanol = int(nrOfMolecules * 228 / 3716) #based on 14.3% mass ethanol
    nrWater = nrOfMolecules - nrEthanol
    d = boxSize / (nrOfMolecules)**(1/3)
    n = int(nrOfMolecules**(1/3)) + 2
    rand = []
    for i in range(1,n): 
        for j in range(1,n): 
            for k in range(1,n): 
                rand.append([(i+1)*d, (j+1)*d, (k+1)*d])
    global indWater, indEthanol
    indWater = range(0, nrOfMolecules)
    indEthanol = np.random.choice(indWater, size=nrEthanol, replace=False) 
    np.delete(indWater, indEthanol)
    
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"{nrWater*3 + nrEthanol*9}\n")
        outputFile.write("Comment t = 0\n")
        for i in indWater:
            outputFile.write(f"O {rand[i][0]} {rand[i][1]} {rand[i][2]}\n")
            outputFile.write(f"H {rand[i][0]-1} {rand[i][1]} {rand[i][2]}\n")
            outputFile.write(f"H {rand[i][0]} {rand[i][1]-0.3} {rand[i][2]-0.9}\n")
        for i in indEthanol:
            outputFile.write(f'C {rand[i][0] + -1} {rand[i][1] + -1} {rand[i][2] + 0.5} \n')
            outputFile.write(f'H {rand[i][0] + -1} {rand[i][1] + -2} {rand[i][2] + 0.2} \n')
            outputFile.write(f'H {rand[i][0] + -2} {rand[i][1] + -1} {rand[i][2] + 1.1} \n')
            outputFile.write(f'H {rand[i][0] + -1} {rand[i][1] + 0} {rand[i][2] + 0.1} \n')
            outputFile.write(f'C {rand[i][0] + 0.5} {rand[i][1] + -1} {rand[i][2] + 0.9} \n')
            outputFile.write(f'H {rand[i][0] + 0.5} {rand[i][1] + -0.2} {rand[i][2] + 1.7} \n')
            outputFile.write(f'H {rand[i][0] + 0.5} {rand[i][1] + -1.7} {rand[i][2] + 1.7} \n')
            outputFile.write(f'O {rand[i][0] + 2} {rand[i][1] + -1} {rand[i][2] + 0.4} \n')
            outputFile.write(f'H {rand[i][0] + 2.5} {rand[i][1] + -0.5} {rand[i][2] + 0.2} \n')


def writeConfig(type, boxSize):
    ''' Writes topology & initial xyz for given system
    
    Number of molecules is based on densities at T = 298.15K '''
    boxSize /= 10 #given in angstrom, computed here with nm
    if (type == 'water'): 
        nrOfMolecules = int(33.328*(boxSize**3))
        writeWaterTopology(nrOfMolecules)
        writeWaterXYZ(nrOfMolecules, boxSize)
    elif (type == 'ethanol'):
        nrOfMolecules = int(10.272*boxSize**3)
        writeEthanolTopology(nrOfMolecules)
        writeEthanolXYZ(nrOfMolecules, boxSize) 
    elif(type == 'mixture'):
        nrOfMolecules = int(29.728*boxSize**3)
        writeMixtureTopology(nrOfMolecules)
        writeMixtureXYZ(nrOfMolecules, boxSize)
    else: 
        print('Type unknown, try again')
            
writeConfig('water', 30)
