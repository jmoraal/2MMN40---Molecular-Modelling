# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 16:30:44 2020

@author: s161981
"""

import numpy as np


lines = []
firstColumn = []

with open("Methane.xyz", "r") as inputFile:
    for line in inputFile: 
        splittedLine = line.split()
        firstColumn.append(splittedLine[0])
        lines.append(splittedLine[1:4])
    
nrOfAtoms = int(firstColumn[0])
timeSteps = int(np.floor(len(lines)/ nrOfAtoms))
#This works because nr of lines not describing positions is much smaller than nrOfAtoms


atomPositions = []
atomTypes = []

for i in range(timeSteps):
    atomTypes.append(firstColumn[(2+(2+nrOfAtoms)*i):((2+nrOfAtoms)*(i+1))])
    atomPositions.append(lines[(2+(2+nrOfAtoms)*i):((2+nrOfAtoms)*(i+1))])
    
atomPositions = np.asarray(atomPositions).astype(np.float)
'''
vragen:
    - lijkt me handig alles van een tijdstip in een keer als
      matrix op te slaan (soort slice van de totale positiematrix),
      kan dat of is het niet efficienter dan per regel?
    - We weten niet van tevoren hoeveel timesteps er zijn, 
      maar dat zou wel dingen makkelijker maken. Is het het 
      waard om aan het begin naar de bestandslengte te kijken?
    - Is atomtypes of nrOfAtoms ooit niet constant?
'''