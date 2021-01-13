#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 11:32:07 2021

@author: Ingrid
"""

"""
DEPENDS ON: vpython
Install (normal): pip3 install vpython
Install (anaconda): conda install -c vpython vpython
Website: https://www.glowscript.org/docs/VPythonDocs/index.html """


from vpython import *
scene = canvas(width=1200, height=800)
atom_radius = 0.3
scene.center = vector(0,0,0)
axes = [vector(1,0,0), vector(0,1,0), vector(0,0,1)]
scene.caption= """
To rotate "camera", drag with right button or Ctrl-drag.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
  On a two-button mouse, middle is left + right.
To pan left/right and up/down, Shift-drag.
"""
atoms = [vector(0.5,2,0),
                vector(1,1,0),
                vector(2.5,1,0),
                vector(3,0.1,0.5)]

rij = atoms[1]-atoms[0]
rjk = atoms[2]-atoms[1]
rkl = atoms[3]-atoms[2]
ro = 0.5*atoms[1]+0.5*atoms[2]
rok = atoms[2] - ro

n1 = (rij).cross(rjk)
n2 = (rjk).cross(rkl)

diri = n1
dirl = -n2
dirk = (-(rok).cross(dirl) + 0.5*(rkl).cross(dirl) + 0.5*(rij).cross(diri)).cross(rok)
dirj = - diri -dirk -dirl

print(diri + dirj + dirk + dirl)
roi = atoms[0] - ro
roj = atoms[1] - ro
rok = atoms[2] - ro
rol = atoms[3] - ro
print((roi).cross(diri) + (roj).cross(dirj) + (rok).cross(dirk) + (rol).cross(dirl))

for at in atoms:
    atom = sphere()
    atom.pos = at
    atom.radius = atom_radius
    atom.color = vector(0,0.58,0.69)
    
for i in range(3):
    curve(atoms[i], atoms[i+1], radius=0.05)
    
# rt = shapes.rectangle(width=3, height=2.6)
# obj = extrusion(path=[vec(0,0,0), vec(0,0,-0.01)], shape=rt, opacity=0.2, color=color.blue)
# obj.rotate(angle=0.1*pi, axis=vector(1,0,0))
# lengthScaler = 0.8

pointer = arrow(pos=atoms[0], axis=diri, shaftwidth=0.05, color=color.red,opacity=0.5)
pointer = arrow(pos=atoms[1], axis=dirj, shaftwidth=0.05, color=color.red,opacity=0.5)
pointer = arrow(pos=atoms[2], axis=dirk, shaftwidth=0.05, color=color.red,opacity=0.5)
pointer = arrow(pos=atoms[3], axis=dirl, shaftwidth=0.05, color=color.red,opacity=0.5)