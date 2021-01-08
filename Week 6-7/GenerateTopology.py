# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:54:00 2021

@author: s161981
"""

### Water ###
def writeWaterTopology(nrOfMolecules, outputFileName = ""): 
    if (len(outputFileName) == 0): 
        outputFileName = 'Water' + str(nrOfMolecules) + 'Topology.txt'
    
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
def writeEthanolTopology(nrOfMolecules, clearFile = True, outputFileName = ""): 
    if (len(outputFileName) == 0): 
        outputFileName = 'Ethanol' + str(nrOfMolecules) + 'Topology.txt'
    
    if clearFile: 
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
def writeMixtureTopology(nrOfMolecules): 
    outputFileName = 'mixture' + str(nrOfMolecules) + 'Topology.txt'
    nrEthanol = int(nrOfMolecules * 228 / 3716) #based on 14.3% mass ethanol
    nrWater = nrOfMolecules - nrEthanol
    
    #This of course doesn' work, as it first creates the entire water topology, then the entire ethanol.
    # writeWaterTopology(nrWater, outputFileName = outputFileName)
    # writeEthanolTopology(nrEthanol, clearFile = False, outputFileName = outputFileName)
    
    #Ugly, but works:
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


writeWaterTopology(150)
writeEthanolTopology(50)
writeMixtureTopology(200)