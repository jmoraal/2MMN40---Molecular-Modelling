# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:37:52 2020

@author: s161981
"""


import numpy as np


lines = []

with open("HydrogenSmall.xyz", "r") as inputFile:
    for line in inputFile: 
        lines.append(line)
    
nrOfAtoms = int(lines[0])
timeSteps = int(np.floor(len(lines)/ nrOfAtoms))
#This works because nr of lines not describing positions is much smaller than nrOfAtoms


atomPositions = np.empty([nrOfAtoms,3,timeSteps])

for i in range(timeSteps):
    for j in range(nrOfAtoms):
        temp = lines[i * nrOfAtoms + j + 2]
        splittedLine = temp.split()
        atomPositions[j,:,i] = float(splittedLine)

#atomPositions[:,:,i] = lines[(2+nrOfAtoms*i):(nrOfAtoms*(i+1)+2)]#[1:4]    
'''    
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