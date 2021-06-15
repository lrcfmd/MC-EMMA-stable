from emma_stable.mc_emma import *
import numpy
""" hexagonal structures """

H_rep=[1,1,1]
H_ap=5.8
Hc=2.45

#Make A blocks
H_A1=Atoms("BaO3",pbc=[1,1,1])
H_A1.set_cell(([H_ap,0.,0.,],[-H_ap/2,H_ap/1.154814,0.],[0.,0.,Hc]))
H_A1.set_scaled_positions([[0.,0.,0.5],[0.,0.5,0.5],[0.5,0.,0.5],[0.5,0.5,0.5]])
H_A1=H_A1.repeat(H_rep)
#view(H_A1)

H_A2=Atoms("BaO",pbc=[1,1,1])
H_A2.set_cell(([H_ap,0.,0.,],[-H_ap/2,H_ap/1.154814,0.],[0.,0.,Hc+0.6]))
H_A2.set_scaled_positions([[0.,0.,0.5],[2./3.,1./3.,0.5]])
H_A2=H_A2.repeat(H_rep)
#view(H_A2)

H_A3=Atoms("BaO",pbc=[1,1,1])
H_A3.set_cell(([H_ap,0.,0.,],[-H_ap/2,H_ap/1.154814,0.],[0.,0.,Hc+0.6]))
H_A3.set_scaled_positions([[0.,0.,0.5],[1./3.,2./3.,0.5]])
H_A3=H_A3.repeat(H_rep)
#view(H_A3)

#Tm site atoms
#Set 1
H_Tm1a=Atoms("In",pbc=[1,1,1])
H_Tm1a.set_cell(([H_ap,0.,0.],[-H_ap/2,(H_ap/2)*1.75,0.],[0.,0.,0.0001]))
H_Tm1a.set_scaled_positions([[0.000,0.000,0.5]])
H_Tm1a=H_Tm1a.repeat(H_rep)
#view(H_Tm1a)

H_Tm1b=Atoms("Al",pbc=[1,1,1])
H_Tm1b.set_cell(([H_ap,0.,0.,],[-H_ap/2,(H_ap/2)*1.75,0.],[0.,0.,0.0001]))
H_Tm1b.set_scaled_positions([[0.000,0.000,0.5]])
H_Tm1b=H_Tm1b.repeat(H_rep)
#view(H_Tm1b)

AMods=[H_A1,H_A2,H_A3]
BMods=[H_Tm1a,H_Tm1b]


mc(
startup=3,
sl=[8,16],
A_Mods={'hexagonal':AMods},
B_Mods={'hexagonal':BMods},
output='full',
charges={"O":-2.,"Ba":2.,"In":3.,"Al":3.},
kwds=['opti conv','opti conp'],
opts=[
	['\nlibrary library.lib\ndump temp.res\ntime 15 minuets\nmaxcyc 100\nstepmax 0.01'],
	['\nlibrary library.lib\ndump temp.res\ntime 15 minutes\nmaxcyc 5000\nlbfgs_order 5000'],
	],
shel=["Ba","O"],
lib = 'library.lib',
smax=1,
rmax=96,
smin=500,
ratt_dist=0.1,
composition={"Ba":2.,"In":1.,"Al":1.,"O":5.},
red_T=0.01,
delay = 0,
check='T',
graph_out='T',
pert=(
['T1',1],
['T2',1],
['T3',1],
['T4',1],
['T5',1],
['T6',1],
),
)



