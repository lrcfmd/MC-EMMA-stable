from emma_stable.mc_emma import *
import numpy
import sys
""" Cubic structures """
rep=[1,1,1]
ap=8 
A1=read("A1.cif")
A2=read("A2.cif")
A3=read("A3.cif")
A4=read("A4.cif")
A5=read("A5.cif")
A6=read("A6.cif")
A7=read("A7.cif")
A8=read("A8.cif")
A9=read("A9.cif")

B1=read("B1.cif")
B2=read("B2.cif")
B3=read("B3.cif")
B4=read("B4.cif")
B5=read("B5.cif")
B6=read("B6.cif")
B7=read("B7.cif")
B8=read("B8.cif")
B9=read("B9.cif")

AMods=[A1,A2,A3,A4,A5,A6,A7,A8,A9]
BMods=[B1,B2,B3,B4,B5,B6,B7,B8,B9]

mc(
#gulp parameters
kwds=['opti conv c6','opti conp c6'],
opts=[
	['\nlibrary lib2.lib\ndump temp.res\ntime 15 minuets\nmaxcyc 100\nstepmax 0.01'],
	['\nlibrary lib2.lib\ndump temp.res\ntime 15 minutes\nmaxcyc 5000\nlbfgs_order 5000'],
	],
shel=[""],
lib = 'lib2.lib',

#compositionparameters
A_Mods={'orthorhombic':AMods},
B_Mods={'orthorhombic':BMods},
charges={"O":-2,"Sr":2,"Ti":4},
composition={"Sr":1.0,"Ti":1.0,"O":3.0},
sl=[1,2,3,4],

#mc run parameters
startup=2,
smax=7,
rmax=10,
smin=100,
ratt_dist=0.2,
red_T=0.01,
delay = 0,
graph_out='T',
pert=(
['T1',29],
['T2',21],
['T3',21],
['T4',14],
['T5',7],
['T6',7],),
kmax=2000000,
tmax=200000,

)

