# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981
"""

import numpy as np

timeSteps = 0

with open("HydrogenSmall.xyz", "r") as inputFile:
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
vragen:
    - lijkt me handig alles van een tijdstip in een keer als
      matrix op te slaan (soort slice van de totale positiematrix),
      kan dat of is het niet efficienter dan per regel?
    - We weten niet van tevoren hoeveel timesteps er zijn, 
      maar dat zou wel dingen makkelijker maken. Is het het 
      waard om aan het begin naar de bestandslengte te kijken?

'''