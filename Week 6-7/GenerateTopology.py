# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:54:00 2021

@author: s161981
"""

### Water ###
nrOfMolecules = 150
outputFileName = 'Water150Topology.txt'

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
# nrOfMolecules = 450
# outputFileName = 'Water450Topology.txt'

# with open(outputFileName, "w") as outputFile: # clear file
#     outputFile.write("") 
# with open(outputFileName, "a") as outputFile:
#     outputFile.write(f"molecules {nrOfMolecules}\n")
#     for i in range(0,nrOfMolecules):
#         outputFile.write(f"{9*i} {9*i+1} {9*i+1} {9*i+1} {9*i+1} {9*i+1} {9*i+1} {9*i+1} {9*i+1}\n")
    
#     outputFile.write(f"bonds {nrOfMolecules*2}\n")
#     for i in range(0,nrOfMolecules):
#         outputFile.write(f"{3*i} {3*i+1} 5024.16 0.9572 \n")
#         outputFile.write(f"{3*i} {3*i+2} 5024.16 0.9572 \n")
    
#     outputFile.write(f"angles {nrOfMolecules*2}\n")
#     for i in range(0,nrOfMolecules):
#         outputFile.write(f"{3*i+1} {3*i} {3*i+2} 628.02 1.8242181 \n")
        
#     outputFile.write(f"LJ {nrOfMolecules*3}\n")
#     for i in range(0,nrOfMolecules):
#         outputFile.write("3.15061 0.66386\n")
#         outputFile.write("0 0\n")
#         outputFile.write("0 0\n") 