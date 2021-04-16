from ase import *
import pickle
import os
import sys
import numpy
from ase.io import *

def write_structures():
    if not os.path.isdir("structures"):
        os.mkdir("structures")
    all_structures=pickle.load(open("all_structures.p",'rb'))
    keys=list(all_structures.keys())
    for i in range(len(keys)):
        label=str("structures/"+str("{0:06d}".format(keys[i]))+".cif")
        write(label,all_structures[keys[i]])   		
write_structures()
