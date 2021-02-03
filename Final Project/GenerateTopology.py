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
            # outputFile.write(f"{9*i +1} {9*i + 0} 2845.12 1.09 \n" )
            # outputFile.write(f"{9*i +2} {9*i + 0} 2845.12 1.09 \n" )
            # outputFile.write(f"{9*i +3} {9*i + 0} 2845.12 1.09 \n" )
            # outputFile.write(f"{9*i +0} {9*i + 4} 2242.624 1.529 \n" )
            # outputFile.write(f"{9*i +5} {9*i + 4} 2845.12 1.09 \n" )
            # outputFile.write(f"{9*i +6} {9*i + 4} 2845.12 1.09 \n" )
            # outputFile.write(f"{9*i +4} {9*i + 7} 2677.76 1.41 \n" )
            # outputFile.write(f"{9*i +7} {9*i + 8} 4627.5 0.945 \n" )
            
            outputFile.write(f"{9*i +0} {9*i + 4} 2242.624 1.529 \n" )
            outputFile.write(f"{9*i +0} {9*i + 2} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +0} {9*i + 3} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +0} {9*i + 1} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +5} {9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +6} {9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{9*i +4} {9*i + 7} 2677.76 1.41 \n" )
            outputFile.write(f"{9*i +7} {9*i + 8} 4627.5 0.945 \n" )
        
        outputFile.write(f"angles {nrOfMolecules*13}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f'{9*i +1} {9*i + 0} {9*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{9*i +2} {9*i + 0} {9*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{9*i +3} {9*i + 0} {9*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{9*i +3} {9*i + 0} {9*i + 2} 276.144 1.88146 \n' )
            outputFile.write(f'{9*i +3} {9*i + 0} {9*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{9*i +2} {9*i + 0} {9*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{9*i +5} {9*i + 4} {9*i + 6} 276.144 1.88146 \n' )
            outputFile.write(f'{9*i +0} {9*i + 4} {9*i + 6} 313.8 1.93208 \n' )
            outputFile.write(f'{9*i +0} {9*i + 4} {9*i + 5} 313.8 1.93208 \n' )
            outputFile.write(f'{9*i +0} {9*i + 4} {9*i + 7} 414.4 1.91114 \n' )
            outputFile.write(f'{9*i +4} {9*i + 7} {9*i + 8} 460.24 1.89368 \n' )
            outputFile.write(f'{9*i +5} {9*i + 4} {9*i + 7} 292.88 1.91114 \n' )
            outputFile.write(f'{9*i +6} {9*i + 4} {9*i + 7} 292.88 1.91114 \n' )

    
        outputFile.write(f"dihedrals {nrOfMolecules*12}\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f'{9*i +1} {9*i + 0} {9*i + 4} {9*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{9*i +2} {9*i + 0} {9*i + 4} {9*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{9*i +3} {9*i + 0} {9*i + 4} {9*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{9*i +1} {9*i + 0} {9*i + 4} {9*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{9*i +2} {9*i + 0} {9*i + 4} {9*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{9*i +3} {9*i + 0} {9*i + 4} {9*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{9*i +1} {9*i + 0} {9*i + 4} {9*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{9*i +2} {9*i + 0} {9*i + 4} {9*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{9*i +3} {9*i + 0} {9*i + 4} {9*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{9*i +0} {9*i + 4} {9*i + 7} {9*i + 8} -0.4431 3.83255 0.72801 -4.11705 \n' )
            outputFile.write(f'{9*i +5} {9*i + 4} {9*i + 7} {9*i + 8} 0.9414 2.8242 0 -3.7656 \n' )
            outputFile.write(f'{9*i +6} {9*i + 4} {9*i + 7} {9*i + 8} 0.9414 2.8242 0 -3.7656 \n' )

        
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
                outputFile.write(f"{nrWater * 3 + 9*i + j} ")
            outputFile.write("\n")
        
        outputFile.write(f"bonds {nrWater*2 + nrEthanol * 8}\n")
        for i in range(0,nrWater):
            outputFile.write(f"{3*i} {3*i+1} 5024.16 0.9572 \n")
            outputFile.write(f"{3*i} {3*i+2} 5024.16 0.9572 \n")
        for i in range(0,nrEthanol):
            # outputFile.write(f"{nrWater*3 + 9*i +1} {nrWater*3 + 9*i + 0} 2845.12 1.09 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +2} {nrWater*3 + 9*i + 0} 2845.12 1.09 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} 2845.12 1.09 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 4} 2242.624 1.529 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +5} {nrWater*3 + 9*i + 4} 2845.12 1.09 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +6} {nrWater*3 + 9*i + 4} 2845.12 1.09 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +4} {nrWater*3 + 9*i + 7} 2677.76 1.41 \n" )
            # outputFile.write(f"{nrWater*3 + 9*i +7} {nrWater*3 + 9*i + 8} 4627.5 0.945 \n" )
            
            outputFile.write(f"{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 4} 2242.624 1.529 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 2} 2845.12 1.09 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 3} 2845.12 1.09 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 1} 2845.12 1.09 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +5} {nrWater*3 + 9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +6} {nrWater*3 + 9*i + 4} 2845.12 1.09 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +4} {nrWater*3 + 9*i + 7} 2677.76 1.41 \n" )
            outputFile.write(f"{nrWater*3 + 9*i +7} {nrWater*3 + 9*i + 8} 4627.5 0.945 \n" )
        
        outputFile.write(f"angles {nrWater + nrEthanol * 13}\n")
        for i in range(0,nrWater):
            outputFile.write(f"{3*i+1} {3*i} {3*i+2} 628.02 1.8242181 \n")
        for i in range(0,nrEthanol):
            outputFile.write(f'{nrWater*3 + 9*i +1} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +2} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} 292.88 1.89368 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 2} 276.144 1.88146 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +2} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 1} 276.144 1.88146 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +5} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 6} 276.144 1.88146 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 6} 313.8 1.93208 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 5} 313.8 1.93208 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} 414.4 1.91114 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +4} {nrWater*3 + 9*i + 7} {nrWater*3 + 9*i + 8} 460.24 1.89368 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +5} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} 292.88 1.91114 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +6} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} 292.88 1.91114 \n' )
            
        outputFile.write(f"dihedrals {nrEthanol * 12}\n")
        for i in range(0,nrEthanol):
            outputFile.write(f'{nrWater*3 + 9*i +1} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +2} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 5} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +1} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +2} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 6} 0.6276 1.8828 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +1} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +2} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +3} {nrWater*3 + 9*i + 0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} 0.97905 2.93716 0 -3.91622 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +0} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} {nrWater*3 + 9*i + 8} -0.4431 3.83255 0.72801 -4.11705 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +5} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} {nrWater*3 + 9*i + 8} 0.9414 2.8242 0 -3.7656 \n' )
            outputFile.write(f'{nrWater*3 + 9*i +6} {nrWater*3 + 9*i + 4} {nrWater*3 + 9*i + 7} {nrWater*3 + 9*i + 8} 0.9414 2.8242 0 -3.7656 \n' )

        
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

def createGrid(nrOfMolecules, boxSize):
    n = int(nrOfMolecules**(1/3)) + 1
    d = boxSize / n
    rand = []
    for i in range(0,n): 
        for j in range(0,n): 
            for k in range(0,n): 
                rand.append([(i+1/2)*d, (j+1/2)*d, (k+1/2)*d]) #1/2 so that initial distance through box edge is equal to distance within
    return rand


def writeWaterXYZ(nrOfMolecules, boxSize):
    outputFileName = 'Water' + str(boxSize) + 'Initial.xyz'
    global rand
    rand = createGrid(nrOfMolecules, boxSize)
    
    
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"{nrOfMolecules*3}\n")
        outputFile.write("Comment t = 0\n")
        for i in range(0,nrOfMolecules):
            outputFile.write(f"O {rand[i][0]} {rand[i][1]} {rand[i][2]}\n")
            outputFile.write(f"H {rand[i][0]-0.98505} {rand[i][1]+0.001286} {rand[i][2]+0.0099}\n")
            outputFile.write(f"H {rand[i][0]+0.011297} {rand[i][1]-0.32759} {rand[i][2]-0.93691}\n")
            
def writeEthanolXYZ(nrOfMolecules, boxSize):
    outputFileName = 'Ethanol' + str(boxSize) + 'Initial.xyz'
    global rand
    rand = createGrid(nrOfMolecules, boxSize)
    
    with open(outputFileName, "w") as outputFile: # clear file
        outputFile.write("") 
    with open(outputFileName, "a") as outputFile:
        outputFile.write(f"{nrOfMolecules*9}\n")
        outputFile.write("Comment t = 0\n")
        for i in range(0,nrOfMolecules):
            # outputFile.write(f'C {rand[i][0] + -0.846410161513776} {rand[i][1] + -1} {rand[i][2] + 0.666025403784439} \n')
            # outputFile.write(f'H {rand[i][0] + -1.10621778264911} {rand[i][1] + -2} {rand[i][2] + 0.516025403784439} \n')
            # outputFile.write(f'H {rand[i][0] + -0.826794919243112} {rand[i][1] + -1} {rand[i][2] + 1.83205080756888} \n')
            # outputFile.write(f'H {rand[i][0] + -1.19282032302755} {rand[i][1] + 0} {rand[i][2] + 0.466025403784438} \n')
            # outputFile.write(f'C {rand[i][0] + 0.25} {rand[i][1] + -1} {rand[i][2] + -0.433012701892219} \n')
            # outputFile.write(f'H {rand[i][0] + 0.942820323027551} {rand[i][1] + -0.2} {rand[i][2] + -0.0330127018922192} \n')
            # outputFile.write(f'H {rand[i][0] + 0.942820323027551} {rand[i][1] + -1.7} {rand[i][2] + -0.0330127018922192} \n')
            # outputFile.write(f'O {rand[i][0] + 0.566987298107781} {rand[i][1] + -1} {rand[i][2] + -1.98205080756888} \n')
            # outputFile.write(f'H {rand[i][0] + 0.643782217350893} {rand[i][1] + -0.5} {rand[i][2] + -2.5150635094611} \n')


            outputFile.write(f'C {rand[i][0] + 3.326028} {rand[i][1] + 2.09962} {rand[i][2] + 1.673972} \n')
            outputFile.write(f'H {rand[i][0] + 3.92445} {rand[i][1] + 1.46277} {rand[i][2] + 2.328371} \n')
            outputFile.write(f'H {rand[i][0] + 3.99617} {rand[i][1] + 2.6448} {rand[i][2] + 1.00383} \n')
            outputFile.write(f'H {rand[i][0] + 2.671629} {rand[i][1] + 1.46277} {rand[i][2] + 1.07555} \n')
            outputFile.write(f'C {rand[i][0] + 2.5} {rand[i][1] + 3.05946} {rand[i][2] + 2.5} \n')
            outputFile.write(f'H {rand[i][0] + 1.903} {rand[i][1] + 3.70751} {rand[i][2] + 1.842751} \n')
            outputFile.write(f'H {rand[i][0] + 3.157249} {rand[i][1] + 3.70751} {rand[i][2] + 3.09706} \n')
            outputFile.write(f'O {rand[i][0] + 1.658486} {rand[i][1] + 2.27233} {rand[i][2] + 3.341514} \n')
            outputFile.write(f'H {rand[i][0] + 1.12353} {rand[i][1] + 2.88153} {rand[i][2] + 3.8764700} \n')



def writeMixtureXYZ(nrOfMolecules, boxSize):
    global nrEthanol, indEthanol, indWater
    outputFileName = 'Mixture' + str(boxSize) + 'Initial.xyz'
    nrEthanol = int(nrOfMolecules * 228 / 3716) #based on 14.3% mass ethanol
    nrWater = nrOfMolecules - nrEthanol
    rand = createGrid(nrOfMolecules, boxSize)
    indWater = np.array(range(0, nrOfMolecules))
    # indEthanol = np.random.choice(indWater, size=nrEthanol, replace=False)
    print(nrOfMolecules)
    print(nrEthanol)
    print(nrOfMolecules/nrEthanol)
    indEthanol = np.arange(0,nrOfMolecules,round(nrOfMolecules/nrEthanol))[:nrEthanol] # TODO for large mixture SUBSTRACT 2 from round(nrOfMolecules/nrEthanol)-2
    indWater = np.delete(indWater, indEthanol, axis=0)
    
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
            # outputFile.write(f'C {rand[i][0] + -0.95} {rand[i][1] + 0.95} {rand[i][2] + 0.777817459305202} \n')
            # outputFile.write(f'H {rand[i][0] + -1.80710678118655} {rand[i][1] + 0.392893218813452} {rand[i][2] + 0.565685424949238} \n')
            # outputFile.write(f'H {rand[i][0] + -1.15} {rand[i][1] + 1.15} {rand[i][2] + 1.90918830920368} \n')
            # outputFile.write(f'H {rand[i][0] + -0.442893218813453} {rand[i][1] + 1.85710678118655} {rand[i][2] + 0.494974746830583} \n')
            # outputFile.write(f'C {rand[i][0] + 0} {rand[i][1] + 0} {rand[i][2] + 0} \n')
            # outputFile.write(f'H {rand[i][0] + 0.965685424949238} {rand[i][1] + 0.165685424949238} {rand[i][2] + 0.565685424949238} \n')
            # outputFile.write(f'H {rand[i][0] + -0.0949747468305832} {rand[i][1] + -0.894974746830583} {rand[i][2] + 0.565685424949238} \n')
            # outputFile.write(f'O {rand[i][0] + 0.5} {rand[i][1] + -0.5} {rand[i][2] + -1.41421356237309} \n')
            # outputFile.write(f'H {rand[i][0] + 1.00355339059327} {rand[i][1] + -0.296446609406726} {rand[i][2] + -1.90918830920368} \n')
            
            
            xshift = -3
            yshift = -3
            zshift = -2
            scale = 0.9
            outputFile.write(f'C {rand[i][0] + scale*(xshift + 3.326028)} {rand[i][1] + scale*(yshift + 2.09962)} {rand[i][2] + scale*(zshift + 1.673972)} \n')
            outputFile.write(f'H {rand[i][0] + scale*(xshift + 3.92445)} {rand[i][1] + scale*(yshift + 1.46277)} {rand[i][2] + scale*(zshift + 2.328371)} \n')
            outputFile.write(f'H {rand[i][0] + scale*(xshift + 3.99617)} {rand[i][1] + scale*(yshift + 2.6448)} {rand[i][2] + scale*(zshift + 1.00383)} \n')
            outputFile.write(f'H {rand[i][0] + scale*(xshift + 2.671629)} {rand[i][1] + scale*(yshift + 1.46277)} {rand[i][2] + scale*(zshift + 1.07555)} \n')
            outputFile.write(f'C {rand[i][0] + scale*(xshift + 2.5)} {rand[i][1] + scale*(yshift + 3.05946)} {rand[i][2] + scale*(zshift + 2.5)} \n')
            outputFile.write(f'H {rand[i][0] + scale*(xshift + 1.903)} {rand[i][1] + scale*(yshift + 3.70751)} {rand[i][2] + scale*(zshift + 1.842751)} \n')
            outputFile.write(f'H {rand[i][0] + scale*(xshift + 3.157249)} {rand[i][1] + scale*(yshift + 3.70751)} {rand[i][2] + scale*(zshift + 3.09706)} \n')
            outputFile.write(f'O {rand[i][0] + scale*(xshift + 1.658486)} {rand[i][1] + scale*(yshift + 2.27233)} {rand[i][2] + scale*(zshift + 3.341514)} \n')
            outputFile.write(f'H {rand[i][0] + scale*(xshift + 1.12353)} {rand[i][1] + scale*(yshift + 2.88153)} {rand[i][2] + scale*(zshift + 3.8764700)} \n')


def writeConfig(type, boxSize):
    ''' Writes topology & initial xyz for given box size (in Angstrom)
    
    Number of molecules is based on densities at T = 298.15K '''
    if (type == 'water'): 
        nrOfMolecules = int(33.328*((boxSize/10)**3))
        writeWaterTopology(nrOfMolecules, boxSize)
        writeWaterXYZ(nrOfMolecules, boxSize)
    elif (type == 'ethanol'):
        nrOfMolecules = int(10.272*((boxSize/10)**3))
        writeEthanolTopology(nrOfMolecules, boxSize)
        writeEthanolXYZ(nrOfMolecules, boxSize) 
    elif(type == 'mixture'):
        nrOfMolecules = int(29.728*((boxSize/10)**3))
        writeMixtureTopology(nrOfMolecules, boxSize)
        writeMixtureXYZ(nrOfMolecules, boxSize)
    else: 
        print('Type unknown, try again')

def findSize(approxSize, substance):
    ''' returns boxsizes so that all grid positions are filled '''
    sizes = np.arange(approxSize - 3.0, approxSize + 3.0, 0.01)
    d = []
    for boxSize in sizes: 
        if (substance == 'water'):
            nrOfMolecules = int(33.328*((boxSize/10)**3)) # water
        elif(substance == 'ethanol'): 
            nrOfMolecules = int(10.272*((boxSize/10)**3)) # ethanol
        elif(substance == 'mixture'):
            nrOfMolecules = int(29.728*((boxSize/10)**3)) # mixture
        else: 
            print('Substance not known')            
        
        n = int(nrOfMolecules**(1/3)) + 1
        if (n**3 - nrOfMolecules ==0): 
            d.append(boxSize)
    return d

#print(findSize(50,'mixture'))

writeConfig('water', 31.08)
# writeConfig('ethanol', 32.22)
# writeConfig('mixture', 51.65)
# writeConfig('mixture', 32.29)

