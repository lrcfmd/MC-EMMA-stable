from emma_stable.all import *
import sys
from random import choice
from random import shuffle
import numpy
import time

def make_random_structure(AMods,BMods,composition,sl,testing=''):
    ## bin modules by their composition
    if testing==True:
    	 pass
    	 #print("\nAMods\n\n",AMods,"\nBMods\n\n",BMods,"\ncomposition\n\n",composition,"\nsl\n\n",sl)

    acomp={}
    bcomp={}
    for i in range(len(BMods)):
        try:
            bcomp[BMods[i].get_chemical_formula()].append(i)
        except KeyError:
            bcomp[BMods[i].get_chemical_formula()]=[i]
    
    for i in range(len(AMods)):
        try:
            acomp[AMods[i].get_chemical_formula()].append(i)
        except KeyError:
            acomp[AMods[i].get_chemical_formula()]=[i]
            
    param1 = list(acomp.keys()) ### A layer compositions
    param2 = list(bcomp.keys()) ### B layer compositions
    
    #print("\n\n",acomp)
    #print("\n\n",bcomp)
    
    if testing == True:
    	 pass
    
    # make sure we've got a list of the unique elements in the calculation
    
    unique_elements=[]
    for i in range(len(AMods)):
    	 atoms=AMods[i]
    	 for j in range(len(atoms)):
    	 	 if atoms[j].symbol not in unique_elements:
    	 	 	 unique_elements.append(atoms[j].symbol)
    
    for i in range(len(BMods)):
    	 atoms=BMods[i]
    	 for j in range(len(atoms)):
    	 	 if atoms[j].symbol not in unique_elements:
    	 	 	 unique_elements.append(atoms[j].symbol)
    
    # collect a module set
    module_set={'A':[],'B':[]} #the two module sequences we'll form into a structure
    #print(composition)
    make = False
    l = choice(sl)
    #l=10
    #print(l)
    for i in range(l):
        mod=choice(param1) # choose A mod to use
    	  #see how many equivalent versions of it there are
        if len(acomp[mod]) == 1:
            module_set['A'].append(acomp[mod][0])
    	 	 
        if len(acomp[mod]) >1:
        	   module_set['A'].append(choice(acomp[mod]))

        mod=choice(param2) # choose B mod to use
        if len(bcomp[mod]) == 1:
            module_set['B'].append(bcomp[mod][0])
    	 	 
        if len(bcomp[mod]) >1:
        	   module_set['B'].append(choice(bcomp[mod]))

    #test to see if the module set has the target composition
    symbols={}
    for i in unique_elements:
    	 symbols[i]=0.
    
    for i in range(l): #go through and add the composition into symbols
        atoms=AMods[module_set['A'][i]].copy()
        syms=atoms.get_chemical_symbols()
        for j in range(len(syms)):
        	  symbols[syms[j]]+=1
        
        atoms=BMods[module_set['B'][i]].copy()
        syms=atoms.get_chemical_symbols()
        for j in range(len(syms)):
        	  symbols[syms[j]]+=1
    
    #print(symbols)
    #check to see if any elements have a value of zero
    keys=list(symbols.keys())
    
    #find element with smallest value
    nums=[]
    for i in range(len(keys)):
        if symbols[keys[i]]==0:
            struct=False
            return struct # attempt failed, return to main code
        nums.append(symbols[keys[i]])

    smallest=keys[nums.index(min(nums))]
    smallest=symbols[smallest]
    #print(smallest)
    new_symbols={}
    
    for i in range(len(keys)):
        new_symbols[keys[i]]=float(symbols[keys[i]]/smallest)
    
    if testing == True:
    	 print(composition)
    	 print(new_symbols)
    	 print(new_symbols == composition)
    	 sys.exit()
    	 #print(module_set)
    	 
    #print(new_symbols)
    #print(composition)
    #time.sleep(2)
    if new_symbols == composition:
        struct=numpy.zeros([l,2],dtype=int)
        shuffle(module_set['A'])
        shuffle(module_set['B'])
        for i in range(len(struct)):
        	  struct[i][0]=module_set['A'][i]
        	  struct[i][1]=module_set['B'][i]
        return struct # if we're going to use this method, convert module_set into a struct object then return
    if new_symbols != composition:
        struct=False
        return struct	 	 
    

