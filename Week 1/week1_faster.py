# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981
"""

import numpy as np


lines = []
firstColumn = []
#NDlines = np.empty([1])

with open("Hydrogen.xyz", "r") as inputFile:
    for line in inputFile: 
        splittedLine = line.split()
        firstColumn.append(splittedLine[0])
        lines.append(splittedLine[1:4])
        #NDlines = np.append(NDlines,np.asarray(line.split()))
    
nrOfAtoms = int(firstColumn[0])
timeSteps = int(np.floor(len(lines)/ nrOfAtoms))
#This works because nr of lines not describing positions is much smaller than nrOfAtoms




#atomPositions = np.empty([nrOfAtoms,timeSteps])
atomPositions = []

for i in range(timeSteps):
    #atomPositions[:,:,i] = extractPosToNDarray(lines[(2+nrOfAtoms*i):(nrOfAtoms*(i+1)+2)])
    atomPositions.append(lines[(2+(2+nrOfAtoms)*i):((2+nrOfAtoms)*(i+1))])
    


     

    
     
'''    


def extractPosToNDarray(linesVar,nrOfAtoms):
    array = np.empty([nrOfAtoms,3])
    for k in range(nrOfAtoms):
        temp = linesVar[k].split()
    
    return array



  nrOfAtoms = int(inputFile.readline()) # now still assuming this is constant throughout the file
  atomType = []
  atomPosition = np.empty([nrOfAtoms,3])
  
  for line in inputFile:
      timeSteps += 1 #can probably be replaced
      for i in range(nrOfAtoms):
          line = inputFile.readline()
          splittedLine = line.split()
          atomType += splittedLine[0]
          temp = splittedLine[1:4]
          atomPosition = np.append(atomPosition,temp)
          try: 
              next(inputFile) #skips the comment line, maybe add later 
          except StopIteration: #cuts out of loop if there are no more lines to read
              break
         
      
      
     
atomPosition = atomPosition.reshape(nrOfAtoms,3,timeSteps)
'''


'''
vragen:
    - lijkt me handig alles van een tijdstip in een keer als
      matrix op te slaan (soort slice van de totale positiematrix),
      kan dat of is het niet efficienter dan per regel?
    - We weten niet van tevoren hoeveel timesteps er zijn, 
      maar dat zou wel dingen makkelijker maken. Is het het 
      waard om aan het begin naar de bestandslengte te kijken?

'''