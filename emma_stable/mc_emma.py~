####### at moment, everything appears to be working with both assembly types, with the exception that when startup = 1 probably needs fixing
import sys
from emma_stable.all import *
import numpy
import math
import os
import sys
from random import choice
import random
import time
from itertools import combinations_with_replacement as cwr
import pickle
from emma_stable.gulp import *
from collections import Counter as ct
import glob
from emma_stable.restart_gulp import *
from emma_stable.make_random_structure import *
import pandas

################################################################################################
import shlex
import re
import numpy
from numpy import arccos, pi, dot
from numpy.linalg import norm

def get_cellpar(atoms):
        cell = atoms.cell
        a = norm(cell[0])
        b = norm(cell[1])
        c = norm(cell[2])
        alpha = arccos(dot(cell[1], cell[2])/(b*c))*180./pi
        beta  = arccos(dot(cell[0], cell[2])/(a*c))*180./pi
        gamma = arccos(dot(cell[0], cell[1])/(a*b))*180./pi
        cell = []
        cell=[a,b,c,alpha,beta,gamma]
        return cell
################################################################################################

#function creates all possible permutations of a given set of modules
def make_permutations(in_array,Mod_limit=''):
    from itertools import permutations
    perms=[] # give perms an empty array
    data=permutations(in_array)
    for array in data:
        tmp=array
    if tmp not in perms:
        if Mod_limit != '':
            if numpy.array(tmp).sum()==Mod_limit: 
                perms.append(array)
        if Mod_limit =='':
            perms.append(array)
    return perms

################################################################################################

def get_distances(new_atoms=''):
    distances=[]
    for i in range(len(new_atoms)):
        for x in range(len(new_atoms)):
            if i != x:
                distances.append(new_atoms.get_distance(i,x,mic=True))
    return distances

################################################################################################

def perovskite_stack(AMods,BMods,SL,struct,ratt_dist=''):
    ### function to convert numerical array structure representation to an atoms object ###
    ### as with cubic EMMA building structure starting from first A layer ###
    atoms=AMods[struct[0,0]]
    tmp_atoms=BMods[struct[0,1]].copy()
    tmp_atoms.translate(atoms.cell[2])
    atoms=atoms+tmp_atoms
    atoms.cell[2]=atoms.cell[2]+BMods[struct[0,1]].cell[2]
    for i in range(1,SL):
        tmp_atoms=AMods[struct[i,0]].copy()
        tmp_atoms.translate(atoms.cell[2])
        atoms=atoms+tmp_atoms
        atoms.cell[2]=atoms.cell[2]+AMods[struct[i,0]].cell[2]
        tmp_atoms=BMods[struct[i,1]].copy()
        tmp_atoms.translate(atoms.cell[2])
        atoms=atoms+tmp_atoms
        atoms.cell[2]=atoms.cell[2]+BMods[struct[i,1]].cell[2]
        del tmp_atoms
    if ratt_dist != '':
        atoms.rattle(ratt_dist)
    return atoms

################################################################################################
#### at the moment, think that this will only work with permutations of the same stacking length ###
def check_struct_cubic(struct,perms,accepted):
    result=0
    length=len(perms)
    periodic_perms=[]
    for i in range(0,len(struct[:,0])):
        tmp_perm=numpy.roll(struct,i,axis=0)
        periodic_perms.append(tmp_perm)
        periodic_perms.append(tmp_perm[::-1])
        del tmp_perm
    periodic_perms=numpy.array(periodic_perms)
    for i in range(len(periodic_perms)):
        for x in range(len(perms)):
            if len(perms[x])==len(periodic_perms[i]):
                tmpA=periodic_perms[i]
                tmpB=perms[x]
                if numpy.array_equal(tmpA,tmpB):
                    result=1
    return result

########################################################################################################################################################################################
################################################################################################
#### at the moment, think that this will only work with permutations of the same stacking length ###


def check_struct_hexagonal(struct,perms,accepted,a_dict='',b_dict='',A_MODS='',B_MODS='',ABC=''):
    result=0
    length=len(perms)
    periodic_perms=[]
    ### rules for checking if hexagonal stacking rules are being applied, ok so this does the A-site stacking, perhaps move this into the structure generation sections?
    if ABC == 'T':
        la=list(a_dict.values())
        lb=list(b_dict.values())

        tmp_perm=struct.copy()
        for i in range(len(tmp_perm)):
            neighbours=[]
            try:
                loc1 = tmp_perm[i][0]
                loc2 = tmp_perm[i+1][0]  
            except(IndexError):
                loc1 = tmp_perm[i][0]
                loc2 = tmp_perm[0][0]
            neighbours.append(loc2)
            loc3=''
            for z in range(len(la)):
                if loc2 in la[z]:
                    loc3=la[z].index(loc2)
                
            for x in range(len(la)):
                if len(la[x]) >> 0:
                    if la[x][loc3] not in neighbours:
                        neighbours.append(la[x][loc3])
    
            if loc1 in neighbours:
                result = 1
                return result


    ### before all of this, need to perform a check to see if the structure obeys hexagonal stacking rules, if ABC == 'T', look to see how it was done in the original EMMA scritps 
    for i in range(0,len(struct[:,0])):
        tmp_perm=numpy.roll(struct,i,axis=0)
        periodic_perms.append(tmp_perm)
        periodic_perms.append(tmp_perm[::-1])

    ### creating first rotation of structures
    rotated_perms=[]

    for i in range(len(periodic_perms)):
        tmp_perm = periodic_perms[i].copy()
        la=list(a_dict.values())
        lb=list(b_dict.values())
        for x in range(len(periodic_perms[i])):
             #### creating the index for the a site module ####
             loc1=''
             loc2=''
             for z in range(len(la)):
                if periodic_perms[i][x][0] in la[z]:
                    loc1=z
                    loc2=la[z].index(int(periodic_perms[i][x][0]))
             #### translating the a-module sequence
             try:
                tmp_perm[x][0]=a_dict[loc1][loc2+1]
             except(IndexError):
                try:
                    tmp_perm[x][0]=a_dict[loc1][loc2-2]
                except(IndexError):### this does not shift the module if the module in question is an empty layer
                    tmp_perm[x][0]=a_dict[loc1][loc2]

             #### creating the index for the b site module ####
             loc1=''
             loc2=''
             for z in range(len(lb)):
                if periodic_perms[i][x][1] in lb[z]:
                    loc1=z
                    loc2=lb[z].index(int(periodic_perms[i][x][1]))
             #### translating the b-module sequence
             try:
                tmp_perm[x][1]=b_dict[loc1][loc2+1]
             except(IndexError):
                try:
                    tmp_perm[x][1]=b_dict[loc1][loc2-2]
                except(IndexError):### this does not shift the module if the module in question is an empty layer
                    tmp_perm[x][1]=b_dict[loc1][loc2]
        rotated_perms.append(tmp_perm)


    #### creating the second translation of each structure #### 
    for i in range(len(periodic_perms)):
        tmp_perm = periodic_perms[i].copy()
        la=list(a_dict.values())
        lb=list(b_dict.values())
        for x in range(len(periodic_perms[i])):
             #### creating the index for the a site module ####
             loc1=''
             loc2=''
             for z in range(len(la)):
                if periodic_perms[i][x][0] in la[z]:
                    loc1=z
                    loc2=la[z].index(int(periodic_perms[i][x][0]))
             #### translating the a-module sequence
             try:
                tmp_perm[x][0]=a_dict[loc1][loc2+2]
             except(IndexError):
                try:
                    tmp_perm[x][0]=a_dict[loc1][loc2-1]
                except(IndexError):### this does not shift the module if the module in question is an empty layer
                    tmp_perm[x][0]=a_dict[loc1][loc2]

             #### creating the index for the b site module ####
             loc1=''
             loc2=''
             for z in range(len(lb)):
                if periodic_perms[i][x][1] in lb[z]:
                    loc1=z
                    loc2=lb[z].index(int(periodic_perms[i][x][1]))
             #### translating the b-module sequence
             try:
                tmp_perm[x][1]=b_dict[loc1][loc2+2]
             except(IndexError):
                try:
                    tmp_perm[x][1]=b_dict[loc1][loc2-1]
                except(IndexError):### this does not shift the module if the module in question is an empty layer
                    tmp_perm[x][1]=b_dict[loc1][loc2]
        rotated_perms.append(tmp_perm)
    for i in range(len(rotated_perms)):
        periodic_perms.append(rotated_perms[i])

    #### in next section, need to create versions in which B and C layers are swapped with each other (moving the B layer with it) ####
    swapped_perms=[]
    for i in range(len(periodic_perms)):
        tmp_perm = periodic_perms[i].copy()
        la=list(a_dict.values())
        lb=list(b_dict.values())
        for x in range(len(periodic_perms[i])):
             #### creating the index for the a site module ####
             loc1=''
             loc2=''
             for z in range(len(la)):###creating the a layer index
                if periodic_perms[i][x][0] in la[z]:
                    loc1=z
                    loc2=la[z].index(int(periodic_perms[i][x][0]))
             loc3=''
             loc4=''
             for z in range(len(lb)):###creating the b layer index
                if periodic_perms[i][x][1] in lb[z]:
                    loc3=z
                    loc4=lb[z].index(int(periodic_perms[i][x][1]))


             if loc2 == 1:
                tmp_perm[x][0]=la[loc1][loc2+1]
             if loc2 == 2:
                tmp_perm[x][0]=la[loc1][loc2-1]
             if loc4 == 1:
                tmp_perm[x][1]=lb[loc3][loc4+1]
             if loc4 == 2:
                 tmp_perm[x][1]=lb[loc3][loc4-1]
        swapped_perms.append(tmp_perm)
    for i in range(len(swapped_perms)):
        periodic_perms.append(swapped_perms[i])    
    periodic_perms=numpy.array(periodic_perms)
    for i in range(len(periodic_perms)):
        for x in range(len(perms)):
            if len(perms[x])==len(periodic_perms[i]):
                tmpA=periodic_perms[i]
                tmpB=perms[x]
                if numpy.array_equal(tmpA,tmpB):
                    result=1
                    break
    return result

########################################################################################################################################################################################

def mc(startup=2,sl='',output='full',A_Mods='',B_Mods='',charges='',ctype='gulp',kwds='',opts='',shel='',ratt_dist='',smax=250,rmax='',lout='lowest_energy_structure.cif',all_out='T',red_T=0.1,check='T',composition='',delay=0,min_bond=1.2,min_bond_num='',
res='T',kmax=10000,smin=1,dump ='T',
start_perm='',start_type='',conj='',conj_steps=0,ABC='T',graph_out='F',
pert=(
['T1',10],
['T2',15],
['T3',15],
['T4',10],
['T5',5],
['T6',2],
)
,r_gulp = ''
,target = ''
,restart = ''
,head = ''
,r_kwrds=''
,r_opts =''
,tmax = 25000
,lib=''
,archive_structures=False
,write_log_file=True
,gulp_command='gulp < gulp.gin > gulp.got'
):

########################################### if r_gulp == True, creating the head.txt file required to restart gulp ########################
    if r_gulp == True:
        z=open(head,'w')
        z.write(str(r_kwrds))
        z.write(str(r_opts))
        z.close()
###note that the pert option needs to be entered in as a tuple of lists, for each list the first element is the pertubation type, the second is it's weighting###
    ###############################################################################################################################################
    ############################################################# removing modules which are not needed ###########################################
    ###############################################################################################################################################
    types=list(A_Mods.keys())
    new_sl={}
    for t in range(len(types)):
        new_sl[types[t]]=[]
        t_i=time.time()
        fsl = []
        nA = len(A_Mods[types[t]])
        nB = len(B_Mods[types[t]])
        AMods=A_Mods[types[t]]
        BMods=B_Mods[types[t]]

    idx=[]
    for i in range(len(AMods)):
        sym=AMods[i].get_chemical_symbols()
        for j in range(len(sym)):
            if not sym[j] in list(composition.keys()):
                if not i in idx:
                    idx.append(i)
    idx.reverse()
    for i in range(len(idx)):
        del AMods[idx[i]]
        
    idx=[]
    for i in range(len(BMods)):
        sym=BMods[i].get_chemical_symbols()
        for j in range(len(sym)):
            if not sym[j] in list(composition.keys()):
                if not i in idx:
                    idx.append(i)
    idx.reverse()
    for i in range(len(idx)):
        del BMods[idx[i]]        
    #for i in range(len(AMods)):
    #    print AMods[i].get_chemical_symbols()        
    #for i in range(len(BMods)):
    #    print BMods[i].get_chemical_symbols()
    #sys.exit()        
    ###############################################################################################################################################    


    if startup != 3:
        all_structures={}
        p_s=0
        allperms={}
        for i in range(len(list(A_Mods.keys()) )):
                allperms[list(A_Mods.keys())[i]]=''
    if dump == 'T':
        if startup ==3:
            o = open("dump.txt",'a')
    if startup == 3:
        acpt=open("accepted_energies.txt",'a')
    if startup != 3:
        if dump =='T':
            o = open("dump.txt",'w')
        acpt=open("accepted_energies.txt",'w')
        if os.path.isdir("steps"):
            old_files = glob.glob("steps/*")
            for i in range(len(old_files)):
                os.remove(old_files[i])
        if archive_structures != True:
            if os.path.isdir("accepted"):
                old_files=glob.glob("accepted/*")
                for i in range(len(old_files)):
                    os.remove(old_files[i])
    
    if startup != 3:
        log_file={'step':[],'starting_modules':[],'modules':[],'move_type':[],'energy':[],
        'mc_outcome':[]}
    	 
    if graph_out =='T' and startup!=3:
        graph=open("energies.txt",'w')
        
    if graph_out =='T' and startup==3:
        graph=open("energies.txt",'a')
    
    t_i = time.time()
    print("\n--------------------------------------------------------------")
    print("|                                                            |")
    print("|              Monte-Carlo EMMA v 4.00                       |")
    print("|                                                            |")
    print("--------------------------------------------------------------") 
    if sl and A_Mods and B_Mods and charges and composition != '' :
        if output == 'full':
            print("input file read")
    else:
        print("ERROR: parameter(s) missing from input file")    
        sys.exit()
    if output == 'full':
        print("requested composition:",composition)
        if archive_structures != True:
            if not os.path.isdir("accepted"):
                os.mkdir("accepted")
    
    if output == 'full':
        print("choosing initial module set")
    if startup ==3:
        if not os.path.isfile("restart.npz"):
            print("restart file not found, starting fresh calculation")
            p_s=0
            allperms={}
            for i in range(len(list(A_Mods.keys()) )):
                    allperms[list(A_Mods.keys())[i]]=''
            startup=2
            
#### counting composisiotns and charges from each layer ####    
    
    if len(types) ==1:
        s_type=types[0]
        start_type=types[0]
    Acharge = {}
    Bcharge = {}
    Acomp_count={}
    Bcomp_count={}
    #### creates translations of hexagonal modules ####
    if 'hexagonal' in types:
        t = types.index('hexagonal')
        tmp_AMods=[]
        a_dict={}
        tmp_BMods=[]
        b_dict={}
        count = 0
        for i in range(len(A_Mods[types[t]])):    
            tmp_AMods.append(A_Mods[types[t]][i])    
            tmp_atoms=A_Mods[types[t]][i].copy()
            tmp_posns=tmp_atoms.get_scaled_positions()
            tmp_cell=get_cellpar(tmp_atoms)
            a_dict[i]=[]
            a_dict[i].append(count)
            count = count+1
            if len(tmp_atoms)>>0:
                tmp_atoms=A_Mods[types[t]][i].copy()
                for x in range(len(tmp_atoms)):
                    tmp_posns[x][0]=tmp_posns[x][0]+1./3.
                    if tmp_posns[x][0] > 1.0:
                        tmp_posns[x][0]=tmp_posns[x][0]-1.
                    tmp_posns[x][1]=tmp_posns[x][1]+2./3.
                    if tmp_posns[x][1] > 1.0:
                        tmp_posns[x][1]=tmp_posns[x][1]-1.
                tmp_atoms.set_scaled_positions(tmp_posns)
                tmp_AMods.append(tmp_atoms)
                a_dict[i].append(count)
                count = count+1
                tmp_atoms=A_Mods[types[t]][i].copy()
                for x in range(len(tmp_atoms)):
                    tmp_posns[x][0]=tmp_posns[x][0]+1./3.
                    if tmp_posns[x][0] > 1.0:
                        tmp_posns[x][0]=tmp_posns[x][0]-1.
                    tmp_posns[x][1]=tmp_posns[x][1]+2./3.
                    if tmp_posns[x][1] > 1.0:
                        tmp_posns[x][1]=tmp_posns[x][1]-1.
                tmp_atoms.set_scaled_positions(tmp_posns)
                tmp_AMods.append(tmp_atoms)
                a_dict[i].append(count)
                count = count +1
        A_Mods[types[t]]=tmp_AMods

        del tmp_atoms
        del tmp_posns
        del tmp_cell
        
        ### making B module translations ### 
        count = 0
        for i in range(len(B_Mods[types[t]])):
            tmp_BMods.append(B_Mods[types[t]][i])    
            tmp_atoms=B_Mods[types[t]][i]
            tmp_posns=tmp_atoms.get_scaled_positions()
            tmp_cell=get_cellpar(tmp_atoms)
            b_dict[i]=[]
            b_dict[i].append(count)
            count = count+1
            if len(tmp_atoms)>>0:
                tmp_atoms=B_Mods[types[t]][i].copy()
                for x in range(len(tmp_atoms)):
                    tmp_posns[x][0]=tmp_posns[x][0]+1./3.
                    if tmp_posns[x][0] > 1.0:
                        tmp_posns[x][0]=tmp_posns[x][0]-1.
                    tmp_posns[x][1]=tmp_posns[x][1]+2./3.
                    if tmp_posns[x][1] > 1.0:
                        tmp_posns[x][1]=tmp_posns[x][1]-1.
                tmp_atoms.set_scaled_positions(tmp_posns)
                tmp_BMods.append(tmp_atoms)

                tmp_atoms=B_Mods[types[t]][i].copy()
                b_dict[i].append(count)
                count = count+1
                for x in range(len(tmp_atoms)):
                    tmp_posns[x][0]=tmp_posns[x][0]+1./3.
                    if tmp_posns[x][0] > 1.0:
                        tmp_posns[x][0]=tmp_posns[x][0]-1.
                    tmp_posns[x][1]=tmp_posns[x][1]+2./3.
                    if tmp_posns[x][1] > 1.0:
                        tmp_posns[x][1]=tmp_posns[x][1]-1.
                tmp_atoms.set_scaled_positions(tmp_posns)
                tmp_BMods.append(tmp_atoms)
                b_dict[i].append(count)
                count = count +1
                
                
        B_Mods[types[t]]=tmp_BMods
        #print "BMods"
        #for i in range(len(B_Mods[types[t]])):
        #        print B_Mods[types[t]][i].get_scaled_positions()
        #        print B_Mods[types[t]][i].get_chemical_symbols()
         #sys.exit()

    for t in range(len(types)):
        new_sl[types[t]]=[]
        t_i=time.time()
        fsl = []
        nA = len(A_Mods[types[t]])
        nB = len(B_Mods[types[t]])
        AMods=A_Mods[types[t]]
        BMods=B_Mods[types[t]]
    
    for t in range(len(types)):
        #count total charge from each layers#
        Acharge[types[t]]=[]
        Bcharge[types[t]]=[]
        Acomp_count[types[t]]={}
        Bcomp_count[types[t]]={}
        nA=len(AMods)
        nB=len(BMods)
        for i in range(nA):
            if len(AMods[i])!=0:
                tmp=AMods[i].get_chemical_symbols()
            tmp_charge=0
            tmp_count=[]
            if len(AMods[i])!=0:
                for x in range(len(tmp)):
                    symbol=(tmp[x])
                    tmp_count.append(symbol)
                    try:
                        tmp_charge=tmp_charge+charges[symbol]
                    except:
                        pass
            Acharge[types[t]].append(tmp_charge)
            Acomp_count[types[t]][i]=tmp_count
        for i in range(nB):
            if len(BMods[i])!=0:
                tmp=BMods[i].get_chemical_symbols()
            tmp_charge=0
            tmp_count=[]
            if len(BMods[i])!=0:
                for x in range(len(tmp)):
                    symbol=(tmp[x])
                    tmp_count.append(symbol)
                    try:
                        tmp_charge=tmp_charge+charges[symbol]
                    except:
                        pass
            Bcharge[types[t]].append(tmp_charge)
            Bcomp_count[types[t]][i]=tmp_count
        amod=list(range(len(AMods)))
        bmod=list(range(len(BMods)))    


###############################################################################################################################################
############################################################# making composition dictionaries##################################################
###############################################################################################################################################

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
    
###############################################################################################################################################
######################## startup from scratch ###################################
#################################################################################
    if rmax == '':
        rmax = max(sl)*max(sl)
        if output == 'full':
            print("no rmax found, set to default:",rmax)

    if startup == 2:
        r=0
        #print "Generating random starting point configuration"
        make = False
        y = 0
        tmax = tmax
        old_A_comp_count = Acomp_count.copy()
        old_B_comp_count = Bcomp_count.copy()
        Acomp_count=old_A_comp_count[types[t]]
        Bcomp_count=old_B_comp_count[types[t]]
        y=0
        while make==False:
            SL=[sl[y]]
            z = 0
            while z < tmax:
                struct=make_random_structure(AMods,BMods,composition,SL)
                if type(struct) != numpy.ndarray:
                    z+=1
                if type(struct) == numpy.ndarray:
                    make = True
                    break
            if z == tmax:
                if y <= len(sl):
                    y+=1
                if y > len(sl):
                	  print("ERROR: failed to generate starting structure")
                	  sys.exit()
                	  
            #print struct    
            #sys.exit()        
#########################################################################################################################    


    #### section in here for hexagonal structures, working out the pools of each translation index for A and B mods
    if types[t]=='hexagonal':
        ### creating pool for A-metal translations, 0 = original, 1= +1/3,+2/3, 2=+2/3,+1/3, 3=vacency layer
        tmp=list(a_dict.keys())
        a_2_dict={}
        a_2_dict[0]=[]
        a_2_dict[1]=[]
        a_2_dict[2]=[]
        a_2_dict[3]=[]        
        for i in range(len(tmp)):
            if len(a_dict[tmp[i]]) >> 1:
                tmp_1=a_dict[tmp[i]]
                a_2_dict[0].append(tmp_1[0])
                a_2_dict[1].append(tmp_1[1])
                a_2_dict[2].append(tmp_1[2])
            if len(a_dict[tmp[i]]) == 1:
                tmp_1=a_dict[tmp[i]]
                a_2_dict[3].append(tmp_1[0])
        ### creating pool for B-metal translations, 0 = original, 1= +1/3,+2/3, 2=+2/3,+1/3, 3=vacency layer
        tmp=list(b_dict.keys())
        b_2_dict={}
        b_2_dict[0]=[]
        b_2_dict[1]=[]
        b_2_dict[2]=[]
        b_2_dict[3]=[]    
        for i in range(len(tmp)):
            if len(b_dict[tmp[i]]) != 1:
                tmp_1=b_dict[tmp[i]]
                b_2_dict[0].append(tmp_1[0])
                b_2_dict[1].append(tmp_1[1])
                b_2_dict[2].append(tmp_1[2])
            if len(b_dict[tmp[i]]) == 1:
                tmp_1=b_dict[tmp[i]]
                b_2_dict[3].append(tmp_1[0])


    #### startup = 1 option, allow the run to start from a user defined array, (also useful for re-starts, i.e. enter in the lowest energy permutation from previous run as starting point for re-start ###
    if startup == 1:
        old_A_comp_count = Acomp_count.copy()
        old_B_comp_count = Bcomp_count.copy()
    	 
        if output == 'full':
            print("starting from user defined configuration")
        if start_perm == '':
            print("ERROR: no input permutation specified, if you want to start from a random permutation, reset startup to 'startup = 2'")
            sys.exit()
        if start_type == '':
            print("ERROR: no stacking type specified")
        
        perms=allperms[start_type]
        struct=numpy.array(start_perm,dtype=int)
        struct_charge=[]
        struct_charge=numpy.array(struct_charge)
        stack=start_type
        AMods=A_Mods[stack][:]
        BMods=B_Mods[stack][:]
        SL = len(struct)
        for i in range(len(struct)):
            struct_charge=numpy.insert(struct_charge,0,Acharge[stack][struct[i,0]])
            struct_charge=numpy.insert(struct_charge,0,Bcharge[stack][struct[i,1]])
        comp_check=0
        gen_comp=[]
        keys=[]
        for i in range(len(struct)):
            tmp=Acomp_count[stack][struct[i,0]]
            for x in range(len(Acomp_count[stack][struct[i,0]])):
                gen_comp.append(Acomp_count[stack][struct[i,0]][x])
            for x in range(len(tmp)):
                if tmp[x] not in keys:
                    keys.append(tmp[x])
                    tmp=Bcomp_count[stack][struct[i,1]]
            for x in range(len(Bcomp_count[stack][struct[i,1]])):
                gen_comp.append(Bcomp_count[stack][struct[i,1]][x])
            for x in range(len(tmp)):
                if tmp[x] not in keys:
                    keys.append(tmp[x])
        trial_comp={}
        for i in range(len(keys)):
            tmpcount=gen_comp.count(keys[i])
            trial_comp[keys[i]]=tmpcount
        min_symb=min(trial_comp,key=trial_comp.get)
        norm = trial_comp[min_symb]
        for sym in range(len(keys)):
            tmp=float(float(trial_comp[keys[sym]])/float(norm))
            trial_comp[keys[sym]]=tmp
        a=trial_comp==composition
        if struct_charge.sum()!=0.0 and a == True:
            print("specified structure not charge neutral, defaulting to random start point")
            startup = 2
            struct = []
        if a == False and struct_charge.sum()==0.0:
            print("specified structure not of specified composition, defaulting to random start point")
            startup = 2
            struct = []
        if a == False and struct_charge.sum()!=0.0:
            print("specified structure not of specified composition and is not charge neutral, defaulting to random start point")
            startup = 2
            struct =[]
        if types[t] == 'hexagonal':
            test =check_struct_hexagonal(struct,perms,accepted=[],a_dict=a_dict,b_dict=b_dict,A_MODS=A_Mods['hexagonal'],B_MODS=B_Mods['hexagonal'],ABC=ABC)
            if test == 1:
                print("specified structure does not have correct hexagonal stacking, defaulting to random starting structure")
                startup = 2
                struct = []
        if types[t] == 'hexagonal':
            comp_check=0
            mods=choice(permitted[SL][0])
            tmpA=[]
            tmpB=[]
            Aset=mods[0]
            Bset=mods[1]
            for y in range(len(Aset)):
                for b in range(Aset[y]):
                    tmpA.append(amod[y])
            zero =0
            one  =0
            two  =0
            three=0
            for y in range(len(tmpA)):
                if tmpA[y] in  a_2_dict[0]:
                    zero = zero +1
                if tmpA[y] in  a_2_dict[1]:
                    one = one +1
                if tmpA[y] in  a_2_dict[2]:
                    two = two +1
                if tmpA[y] in a_2_dict[3]:
                    three=three+1
            test=list([zero,one,two])
            if float(max(test)) > float(len(tmpA)/2):
                comp_check=1
            if float(max(test)) <= float(len(tmpA)/2):
                b = 0
                while b == 0:
                    b = 1
                    numpy.random.shuffle(tmpA)
                    la=list(a_dict.values())
                    lb=list(b_dict.values())
                    tmp_perm=tmpA
                    for i in range(len(tmp_perm)):
                        neighbours=[]
                        try:
                            loc1 = tmp_perm[i]
                            loc2 = tmp_perm[i+1]
                        except(IndexError):
                            loc1 = tmp_perm[i]
                            loc2 = tmp_perm[0]
                        neighbours.append(loc2)
                        loc3=''
                        for z in range(len(la)):
                            if loc2 in la[z]:
                               loc3=la[z].index(loc2)
                            
                        for x in range(len(la)):
                            if len(la[x]) >> 0:
                                if la[x][loc3] not in neighbours:
                                    neighbours.append(la[x][loc3])
                        
                        if loc1 in neighbours:
                            b = 0
            ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
            if comp_check !=1:
                for y in range(len(Bset)):
                    for b in range(Bset[y]):
                        tmpB.append(bmod[y])
                numpy.random.shuffle(tmpB)
                for b in range(len(tmpB)):
                    trans=[0,0]
                    ### identifying the A layer level with the B module
                    if tmpA[b] in b_2_dict[0]:
                        trans[0]=0
                    if tmpA[b] in b_2_dict[1]:
                        trans[0]=1
                    if tmpA[b] in b_2_dict[2]:
                        trans[0]=2
                    if tmpA[b] in b_2_dict[3]:
                        trans[0]=3
                    ### identifying the neihbouring a module
                    if b != max(range(len(tmpB))):
                        if tmpA[b+1] in a_2_dict[0]:
                            trans[1]=0
                        if tmpA[b+1] in a_2_dict[1]:
                            trans[1]=1
                        if tmpA[b+1] in a_2_dict[2]:
                            trans[1]=2
                        if tmpA[b+1] in a_2_dict[3]:
                            trans[1]=3
                    if b+1 == len(tmpB):
                        if tmpA[0] in a_2_dict[0]:
                            trans[1]=0
                        if tmpA[0] in a_2_dict[1]:
                            trans[1]=1
                        if tmpA[0] in a_2_dict[2]:
                            trans[1]=2
                        if tmpA[0] in a_2_dict[3]:
                            trans[1]=3
                    trans=numpy.array(trans)
                    b_mod_idx=0
                    if tmpB[b] in b_2_dict[0]:
                        b_mod_idx=0
                    if tmpB[b] in b_2_dict[1]:
                        b_mod_idx=1
                    if tmpB[b] in b_2_dict[2]:
                        b_mod_idx=2
                    if tmpB[b] in b_2_dict[3]:
                        b_mod_idx=3


                    if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all() == True:
                        trans_1=2
                        if trans_1 == b_mod_idx or b_mod_idx==3:
                            pass
                        if trans_1 != b_mod_idx and b_mod_idx != 3:
                            for x in range(len(list(b_dict.keys()))):
                                if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                    loc1=x
                                    loc2=b_dict[x]
                                    tmpB[b]=b_dict[x][trans_1]
                    if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                        trans_1=1
                        if trans_1 == b_mod_idx or b_mod_idx==3:
                            pass
                        if trans_1 != b_mod_idx and b_mod_idx != 3:
                            for x in range(len(list(b_dict.keys()))):
                                if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                    loc1=x
                                    loc2=b_dict[x]
                                    tmpB[b]=b_dict[x][trans_1]
                    if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                        trans_1=0
                        if trans_1 == b_mod_idx or b_mod_idx==3:
                            pass
                        if trans_1 != b_mod_idx and b_mod_idx != 3:
                            for x in range(len(list(b_dict.keys()))):
                                if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                    loc1=x
                                    loc2=b_dict[x]
                                    tmpB[b]=b_dict[x][trans_1]
            if comp_check !=1:
                for x in range(len(struct)):
                    struct[x,0]=tmpA[x]
                    struct[x,1]=tmpB[x]
                ### ok, now correctly identifying what the b-site translations should be and what they are.... now need to go back and find a way if they are different to change the module for a diffreent translation, but of the same chemical composition... think this should be done using the b_dict dictionary, so you can pull a value out of the same key, but this time with the correct index
    
    #### startup = 2: default option, start from random structure ###
    if startup == 2:
        accept =False
        print("Generating random starting point configuration")
        make = False
        tmax = tmax
        #old_A_comp_count = Acomp_count.copy()
        #old_B_comp_count = Bcomp_count.copy()
        #print Acomp_count
        #print Bcomp_count
        #sys.exit()
        Acomp_count=old_A_comp_count[types[t]]
        Bcomp_count=old_B_comp_count[types[t]]
        n = 0
        while accept == False:
            z = 0
            y = 0
            # initially try to generate a structure, keeping it as small as possible!
            while z <= tmax:
                Acomp_count=old_A_comp_count[types[t]]
                Bcomp_count=old_B_comp_count[types[t]]
                #print "y ",y
                l = sl[y] #make stacking length equal to current trial
                struct=numpy.zeros([l,2],dtype=int)
                for i in range(len(struct)):
                    struct[i][0]=choice(list(range(len(AMods))))
                    struct[i][1]=choice(list(range(len(BMods))))
                    
                struct_charge=[]
                struct_charge=numpy.array(struct_charge)
                SL = len(struct)
                for i in range(len(struct)):
                    struct_charge=numpy.insert(struct_charge,0,Acharge[types[t]][struct[i,0]])
                    struct_charge=numpy.insert(struct_charge,0,Bcharge[types[t]][struct[i,1]])
                
                comp_check=0
                gen_comp=[]
                keys=[]
                for i in range(len(struct)):
                    tmp=Acomp_count[struct[i,0]]
                    for x in range(len(Acomp_count[struct[i,0]])):
                        gen_comp.append(Acomp_count[struct[i,0]][x])
                    for x in range(len(tmp)):
                        if tmp[x] not in keys:
                            keys.append(tmp[x])
                    tmp=Bcomp_count[struct[i,1]]
                    for x in range(len(Bcomp_count[struct[i,1]])):
                        gen_comp.append(Bcomp_count[struct[i,1]][x])
                    for x in range(len(tmp)):
                        if tmp[x] not in keys:
                            keys.append(tmp[x])
                trial_comp={}
                for i in range(len(keys)):
                    tmpcount=gen_comp.count(keys[i])
                    trial_comp[keys[i]]=tmpcount
                min_symb=min(trial_comp,key=trial_comp.get)
                norm = trial_comp[min_symb]
                for sym in range(len(keys)):
                    tmp=float(float(trial_comp[keys[sym]])/float(norm))
                    trial_comp[keys[sym]]=tmp
                a=trial_comp==composition
                
                if a == True:
                    make = True
                    break
                
                z = z+1
                if z == tmax:
                    if sl[y] < max(sl):
                        y = y+1
                        z = 0
                        
                    if sl[y] > max(sl):
                        break
                        
                    if sl[y] == max(sl):
                        if z == tmax:
                            break
            
            if make != True:    
                z = 0
                tmax = 500000
                # Have a second go at maximum length
                while z <= tmax:
                    l = max(sl) #make stacking length equal to current trial
                    struct=numpy.zeros([l,2],dtype=int)
                    for i in range(len(struct)):
                        struct[i][0]=choice(list(range(len(AMods))))
                        struct[i][1]=choice(list(range(len(BMods))))
                        
                    struct_charge=[]
                    struct_charge=numpy.array(struct_charge)
                    SL = len(struct)
                    for i in range(len(struct)):
                        struct_charge=numpy.insert(struct_charge,0,Acharge[types[t]][struct[i,0]])
                        struct_charge=numpy.insert(struct_charge,0,Bcharge[types[t]][struct[i,1]])
                    comp_check=0
                    gen_comp=[]
                    keys=[]
                    for i in range(len(struct)):
                        tmp=Acomp_count[struct[i,0]]
                        for x in range(len(Acomp_count[struct[i,0]])):
                            gen_comp.append(Acomp_count[struct[i,0]][x])
                        for x in range(len(tmp)):
                            if tmp[x] not in keys:
                                keys.append(tmp[x])
                        tmp=Bcomp_count[struct[i,1]]
                        for x in range(len(Bcomp_count[struct[i,1]])):
                            gen_comp.append(Bcomp_count[struct[i,1]][x])
                        for x in range(len(tmp)):
                            if tmp[x] not in keys:
                                keys.append(tmp[x])
                    trial_comp={}
                    for i in range(len(keys)):
                        tmpcount=gen_comp.count(keys[i])
                        trial_comp[keys[i]]=tmpcount
                    min_symb=min(trial_comp,key=trial_comp.get)
                    norm = trial_comp[min_symb]
                    for sym in range(len(keys)):
                        tmp=float(float(trial_comp[keys[sym]])/float(norm))
                        trial_comp[keys[sym]]=tmp
                    a=trial_comp==composition
                    
                    if a == True:
                        make = True
                        break
                    z = z+1
                    if z == tmax:
                        if make == False:
                            print("Error: failed to generate a starting structure")
                            sys.exit()            
            #n = 0
            n_tot=0
            z = 1
            m = 0
            stack = types[t]
            AMods=A_Mods[stack][:]
            BMods=B_Mods[stack][:]
            perms=allperms[stack]
            #generate first structure at random, accecpting the first generated structure which is charge neutral#
            #when the structure type is hexagonal, now need to go through and setup the second array defining the translations
            #struct=numpy.zeros((SL,2),dtype=int)
            #struct_charge=[]
            #struct_charge=numpy.array(struct_charge)
            #print struct
            if types[t] == 'hexagonal':
                comp_check=0
                struct_Amods=[]
                struct_Bmods=[]
                for mod in range(len(struct)):
                    struct_Amods.append(struct[mod][0])
                    struct_Bmods.append(struct[mod][1])
                mods = []
                ModA=()
                ModB=()
                for i in range(len(A_Mods[types[t]])):
                        ModA += (struct_Amods.count(i),)
                for i in range(len(B_Mods[types[t]])):
                        ModB += (struct_Bmods.count(i),)
                mods=[ModA,ModB]
                tmpA=struct_Amods
                tmpB=struct_Bmods
                Aset=mods[0]
                Bset=mods[1]
                #for y in range(len(Aset)):
                #    for b in range(Aset[y]):
                #        tmpA.append(amod[y])
                zero =0
                one  =0
                two  =0
                three=0
                for y in range(len(tmpA)):
                    if tmpA[y] in  a_2_dict[0]:
                        zero = zero +1
                    if tmpA[y] in  a_2_dict[1]:
                        one = one +1
                    if tmpA[y] in  a_2_dict[2]:
                        two = two +1
                    if tmpA[y] in a_2_dict[3]:
                        three=three+1
                test=list([zero,one,two])
                #sys.exit()
                if float(max(test)) > float(len(tmpA)/2):
                    comp_check=1
                if float(max(test)) <= float(len(tmpA)/2):
                    b = 0
                    while b == 0:
                        b = 1
                        numpy.random.shuffle(tmpA)
                        la=list(a_dict.values())
                        lb=list(b_dict.values())
                        tmp_perm=tmpA
                        for i in range(len(tmp_perm)):
                            neighbours=[]
                            try:
                                loc1 = tmp_perm[i]
                                loc2 = tmp_perm[i+1]
                            except(IndexError):
                                loc1 = tmp_perm[i]
                                loc2 = tmp_perm[0]
                            neighbours.append(loc2)
                            loc3=''
                            for z in range(len(la)):
                                if loc2 in la[z]:
                                    loc3=la[z].index(loc2)
                                
                            for x in range(len(la)):
                                if len(la[x]) >> 0:
                                    if la[x][loc3] not in neighbours:
                                        neighbours.append(la[x][loc3])
                               
                            if loc1 in neighbours:
                                b = 0
                #sys.exit()
                ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
                if comp_check !=1:
                    #for y in range(len(Bset)):
                    #    for b in range(Bset[y]):
                    #        tmpB.append(bmod[y])
                    numpy.random.shuffle(tmpB)
                    for b in range(len(tmpB)):
                        trans=[0,0]
                        ### identifying the A layer level with the B module
                        if tmpA[b] in b_2_dict[0]:
                            trans[0]=0
                        if tmpA[b] in b_2_dict[1]:
                            trans[0]=1
                        if tmpA[b] in b_2_dict[2]:
                            trans[0]=2
                        if tmpA[b] in b_2_dict[3]:
                            trans[0]=3
                        ### identifying the neihbouring a module
                        if b != max(range(len(tmpB))):
                            if tmpA[b+1] in a_2_dict[0]:
                                trans[1]=0
                            if tmpA[b+1] in a_2_dict[1]:
                                trans[1]=1
                            if tmpA[b+1] in a_2_dict[2]:
                                trans[1]=2
                            if tmpA[b+1] in a_2_dict[3]:
                                trans[1]=3
                        if b+1 == len(tmpB):
                            if tmpA[0] in a_2_dict[0]:
                                trans[1]=0
                            if tmpA[0] in a_2_dict[1]:
                                trans[1]=1
                            if tmpA[0] in a_2_dict[2]:
                                trans[1]=2
                            if tmpA[0] in a_2_dict[3]:
                                trans[1]=3
                        trans=numpy.array(trans)
                        b_mod_idx=0
                        if tmpB[b] in b_2_dict[0]:
                            b_mod_idx=0
                        if tmpB[b] in b_2_dict[1]:
                            b_mod_idx=1
                        if tmpB[b] in b_2_dict[2]:
                            b_mod_idx=2
                        if tmpB[b] in b_2_dict[3]:
                            b_mod_idx=3
            
            
                        if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all() == True:
                            trans_1=2
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                        if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                            trans_1=1
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                        if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                            trans_1=0
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                if comp_check !=1:
                    for x in range(len(struct)):
                        struct[x,0]=tmpA[x]
                        struct[x,1]=tmpB[x]
                    ### ok, now correctly identifying what the b-site translations should be and what they are.... now need to go back and find a way if they are different to change the module for a diffreent translation, but of the same chemical composition... think this should be done using the b_dict dictionary, so you can pull a value out of the same key, but this time with the correct index
                #print struct    
                #sys.exit()
            
            else:
                tmpA=[]
                tmpB=[]
                struct_modA=()
                struct_modB=()
                for i in range(len(struct)):
                    tmpA.append(struct[i][0])
                    tmpB.append(struct[i][1])
                for i in range(len(A_Mods)):
                    struct_modA+=(tmpA.count(i),)
                for i in range(len(B_Mods)):
                    struct_modB+=(tmpB.count(i),)
                mods=(tmpA,tmpB)
                Aset=mods[0]
                Bset=mods[1]
                #for y in range(len(Aset)):
                #    for b in range(Aset[y]):
                #        tmpA.append(amod[y])
                #for y in range(len(Bset)):
                #    for b in range(Bset[y]):
                #        tmpB.append(bmod[y])
                tmpA=numpy.array(tmpA)
                tmpB=numpy.array(tmpB)
                #numpy.random.shuffle(tmpA)
                #numpy.random.shuffle(tmpB)
                #for i in range(len(struct)):
                #    struct[i,0]=tmpA[i]
                #    struct[i,1]=tmpB[i] 
                struct_charge=[0] #### is this due to me removing the charge counting? as I can't seem to find it now... ####
                struct_charge=numpy.array(struct_charge)
                comp_check=0
            
            #sys.exit()
            ### this is the composition checking section, so this should apply to all structure types ###
            keys=[]
            gen_comp=[]
            Acomp_count=old_A_comp_count.copy()
            Bcomp_count=old_B_comp_count.copy()
            #print Acomp_count
            #print stack
            #sys.exit()
            for i in range(len(struct)):
                tmp=Acomp_count[stack][struct[i,0]]
                for x in range(len(Acomp_count[stack][struct[i,0]])):
                    gen_comp.append(Acomp_count[stack][struct[i,0]][x])
                for x in range(len(tmp)):
                    if tmp[x] not in keys:
                        keys.append(tmp[x])
                tmp=Bcomp_count[stack][struct[i,1]]
                for x in range(len(Bcomp_count[stack][struct[i,1]])):
                    gen_comp.append(Bcomp_count[stack][struct[i,1]][x])
                for x in range(len(tmp)):
                    if tmp[x] not in keys:
                        keys.append(tmp[x])
            trial_comp={}
            for i in range(len(keys)):
                tmpcount=gen_comp.count(keys[i])
                trial_comp[keys[i]]=tmpcount
            min_symb=min(trial_comp,key=trial_comp.get)
            norm = trial_comp[min_symb]
            for sym in range(len(keys)):
                tmp=float(float(trial_comp[keys[sym]])/float(norm))
                trial_comp[keys[sym]]=tmp
                
            a=trial_comp==composition
            if a==False:
                comp_check=1
            
            if comp_check==0:
                test =  0
                if str(s_type) == 'hexagonal':
                    test = check_struct_hexagonal(struct,perms,accepted=[],a_dict=a_dict,b_dict=b_dict,A_MODS=A_Mods['hexagonal'],B_MODS=B_Mods['hexagonal'],ABC=ABC)
                if test == 0:
                    new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
                    if min_bond_num == '':
                        distances=get_distances(new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist))
                        if min(distances) <= min_bond:
                            test = 1
                    #if not min(distances) <= min_bond:
                    if min_bond_num != '':
                        #print "flag2"
### min_bond_nu, should be a parameter to use to specifiy a minimum number of bonds for different species, expressed as a dictionary. e.g. min_bond_num = {"Al":[4,'O',3]}, note that this specifies min of 4 bonds for Al to Oxygen, with a cufoff of 3 Ang.
                        if len(list(min_bond_num.keys())) >> 0:
                            bonds=[]
                            for i in range(len(list(min_bond_num.keys()))):
                                cutoff=float(min_bond_num[list(min_bond_num.keys())[i]][2])
                                num=min_bond_num[list(min_bond_num.keys())[i]][0]
                                s1=list(min_bond_num.keys())[i]
                                #print s1
                                s2=min_bond_num[list(min_bond_num.keys())[i]][1]
                                s1a=[]
                                s2a=[]
                                for a in range(len(new_atoms)):
                                    if new_atoms[a].symbol==s1:
                                        s1a.append(a)
                                    if new_atoms[a].symbol==s2:
                                        s2a.append(a)
                                for x in range(len(s1a)):
                                    rij=[]
                                    for y in range(len(s2a)):
                                        tmp=new_atoms.get_distance(s1a[x],s2a[y],mic=True)
                                        if tmp <= cutoff and tmp >= min_bond:
                                            rij.append(tmp)
                                    bonds.append(rij)
                            for x in range(len(bonds)):
                                #print len(bonds[x])
                                if len(bonds[x]) < num:
                                    test=1
                            
                #print "test",test
                if test == 0:
                    z=0
                    accept = True
                else:
                    z=1
                    make = False
            n = n+1
            #print n
            if n == 50000:
                if len(sl) ==1:
                    print("Error: could not generate charge neutral structure with given input and initial stacking length")
                    sys.exit()
                if len(sl) >>1:
                    while z == 1:
                        m = m+1
                        if m >= 500:
                            print("Error: could not generate charge neutral structure with given input")
                            sys.exit()
                        if SL not in fsl:
                            fsl.append(SL)
                        SL=choice(sl)
                        if SL not in fsl:
                            struct=numpy.zeros((SL,2),dtype=int)
                            n = 0
                            continue
    #sys.exit()
    
    if startup == 3: ## restarting previous run, read in previous: perms, energies, allowed sl & previous lowest energy structure
        try:
        	  log_file=pandas.read_csv('log_file.csv')
        	  log_file=log_file.to_dict(orient='list')
        except:
            log_file={'step':[],'starting_modules':[],'modules':[],'move_type':[],'energy':[],
            'mc_outcome':[]}
        tmp = numpy.load("restart.npz",allow_pickle=True)
        perms=list(tmp['perms'].copy())
        #allperms=tmp_perms.item().copy()
        accepted =tmp['accepted'].copy()
        accepted=list(accepted)
        struct=tmp['lowest'].copy()
        stack=tmp['cstack']
        stack=stack.item()
        #perms=list(allperms[stack])
        #print allperms
        #print "perms",perms
        SL = len(struct)
        #sl=tmp['sl'].copy()
        #sl=sl.item().copy()
        #tmp_s = open("Sl.dat",'r')
        #sl = pickle.load(tmp_s)
        eng = tmp['energies']
        energies = {}
        if archive_structures==True:
            all_structures=pickle.load(open("all_structures.p",'rb'))
        for i in range(len(eng)):
            energies[i]=eng[i,1]
        p_s=max(energies)
        #if res=='T':
        r =int(tmp['r']) 
        tmp = ''
        eng = ''
        old_sl = SL
        tmp=sl
#        for i in range(len(tmp)):
        print("stacking type:",str(types[t]+","), end=' ')
        print("possible stacking lengths=",sl)
        print("restart file read")
        tmp = ''
        old_A_comp_count= Acomp_count.copy()
        old_B_comp_count= Bcomp_count.copy()

    if output == 'full':
        print("Initial stacking is",stack,"with stacking length:",SL)
    atoms = perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
    #view(atoms)
    #sys.exit()
    ### calculation step for initial structure ###
    if ctype == 'gulp':
        if output == 'full':
            print("GULP to be used as calculator")
        write("start_point.cif",atoms)
        s = 0

        if startup !=3:
            atoms,initial_energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=opts,lib=lib,gulp_command=gulp_command)
            if converged == False:
                initial_energy = 1.0E20
                force = 10
                f_type = 'high'
            if graph_out=='T':
                graph.write(str(str(s+p_s).ljust(8)+str('{0:6s}'.format(str("NE = ")))+str('{0:4.8f}'.format(initial_energy/len(atoms)).ljust(35))+str(" eV/atom")))
        #view(atoms)                    
        if startup ==3:
            idx=min(energies, key=energies.get)
            initial_energy=energies[idx]
            idx
############# check in here to make sure that atom - atom distances are longer that min_bond (from input, with default = 1.2 A), to make sure that structures from GULP that were on the verge of collapsing are not accepted #################
    for i in range(len(atoms)):
        tmp = []
        for x in range(len(atoms)):
            tmp.append(x)
        del tmp[i]
        distances=[]
        for x in range(len(tmp)):
            distances.append(atoms.get_distance(i,tmp[x],mic=True))
    if float(min(distances))<= float(min_bond):
        initial_energy = 1.0E20
    if startup != 3:
        try:    
            initial_energy=(initial_energy/len(atoms))    
        except(TypeError,ValueError):
            initial_energy = 1.0e20/len(atom)
    if all_out=='T':
        if archive_structures != True:
            if os.path.isdir("steps") == False:
                os.mkdir("steps")
        if startup != 3:
            #filename = str("steps/step_"+str(s+p_s)+".cif")
            #write(filename,atoms)
            if archive_structures == True:
                all_structures[0]=atoms
            if archive_structures != True:
                write("steps/"+str('{0:0=5n}'.format(s+p_s)+".cif"),atoms)
    t_s=time.time()-t_i
    if startup !=3:
        print("initial configuration energy= %6.4e"%(initial_energy), " eV/atom CPU: %3.1f"%(t_s))
        o.write(str("initial configuration energy= {0: 6.4e} eV/atom, CPU: {1:3.1f}".format(initial_energy,t_s)+"\n\n"))
        if initial_energy < 1.0e+5:
            acpt.write(str("\n"+str(s+p_s).ljust(5)+str('{0:4s}'.format(str("NE = ")))+str('{0:4.8f}'.format(initial_energy).ljust(35))+str(" eV/atom")))
    if startup ==3:
        print("\ncurrent lowest energy = %6.4e"%(initial_energy),"eV/atom, step:",idx,"\n")
        #o.write(str("\ncurrent lowest energy = {0: 6.4e} eV/atom, step: {1:n}".format(initial_energy,idx)+"\n\n"))
    o.flush()
    acpt.flush()
    #if startup != 3:
    #    outfile=open('permitted.dat','wb')
    #    pickle.dump(permitted,outfile,protocol=pickle.HIGHEST_PROTOCOL)


    print("\n*** entering MC loop ***\n")
    accepted=[]
    accepted.append(struct)
    #write("accepted/"+str(s+p_s)+".cif",atoms)
    if archive_structures != True:
        write("accepted/"+str('{0:0=5n}'.format(s+p_s)+".cif"),atoms)
    if startup != 3:
        perms = []
        perms.append(struct)
    l_conf_atoms=atoms
    lowest_energy_structure=atoms.copy()
    lowest_struct=struct.copy()
    lowest_energy=initial_energy
    if startup != 3:
        energies={}    
        energies[(s+p_s)]=initial_energy
    s = 1
    ### generating weighted list to use for the pertubation step, need to extend this section to allow for the weighting to change if the number of rejected structures gets beyond a certain threshold ###
    c_f = []
    for a in range(len(pert)):
        for b in range(pert[a][1]):
            c_f.append(pert[a][0])
    p_log=[]
    #if res != 'T':
    #    r=0
    if res == 'T' and startup !=3:
        r = 0
    k=0
    outcome ='' 
    force = delay    ### force a restriction on the pertubation type, for n trial structures
    f_type = 'low'        ### low  forces swaps to be either type 1 or 2
    ### making the main MC loop ###
    while s <= smax:
        #if outcome=='R':
        struct=accepted[-1]
        if not s+p_s in log_file['step']:
            log_file['step'].append(s+p_s)
            log_file['starting_modules'].append(struct)
        SL=len(struct)
        ### pertubation step to generate new atoms object ###
        if force == 0: ### when enabled types other than 1 & 2, if k reaches threshold, set force = 1, then force p_t to be from types 3 onwards, with equal weighting, or in some cases force to select from only types 1 & 2 only
            p_t=choice(c_f)
        if force != 0 and f_type == 'low':
            p_t=choice(['T1','T2'])
            force = force-1
        if force != 0 and f_type == 'high': ### at the moment, used when r >= 90% rmax, forcing the pertubations to be of the kind the will trial different typs/number of layers
            tmp_f = []
            for i in range(len(pert)):
                #tmp_f = pert[i][0]
                if pert[i][0] in ['T3','T4','T5','T6','T7']:
                    tmp_f.append(pert[i][0])
            if len(tmp_f)>>0:
                p_t = choice(tmp_f)                    
            force = force-1
        if force != 0 and f_type == 'T4':     # what are these here for??
            p_t = 'T4'                    # what are these here for??
            force = force -1              # what are these here for??
        #flag4
        #print p_t


        if p_t == 'T1': ### type 1: two modules should swap place within current structure
            old_sl=SL
            SL = len(struct)
            if types[t]=='hexagonal':
                temp=struct.copy()
                tmpA=temp[:,0]
                #print tmpA
                tmpB=list(set(tmpA))
                accept=False
                if len(tmpB)>=3:
                    trial2=0
                    while accept==False and trial2 <= 100:
                        tmpA=temp[:,0].copy()
                        accept=True ### works on accepting structure as valid, unless proven otehrwise
                        ### first bit, choosing two modules from the A-site sequence and swapping thier positions ###
                        c1=choice(list(range(len(tmpA))))
                        c2=choice(list(range(len(tmpA))))
                        trial = 0
                        while tmpA[c1] == tmpA[c2] and trial <= 50:
                            c2=choice(list(range(len(tmpA))))
                            trial=trial+1
                        #print c1, c2
                        tmpC=[]
                        for i in range(len(tmpA)):
                            tmpC.append(tmpA[i])
                        #print tmpC
                        tmpC[c1]=tmpA[c2].copy()
                        tmpC[c2]=tmpA[c1].copy()
                        #print tmpC
                        ### next bit, working out what the neighbouring modules are, if the neighbours of the the swapped positions are equal to the swaped layers, structure rejected, exception needs to be included here for the eventuality that the modules are vacancy layers... i.e. if len module X = 0.
                        try:
                            n1=tmpC[c1+1]
                        except(IndexError):
                            n1=tmpC[0]
                        n2 = tmpC[c1-1]
                        try:
                            n3=tmpC[c2+1]
                        except(IndexError):
                            n3=tmpC[0]
                        n4 = tmpC[c2-1]
                        if tmpC[c1] == n1 or tmpC[c1] == n2:
                            accept=False
                        if tmpC[c2] == n3 or tmpC[c2] == n4:
                            accept=False
                        ### next bit is checking the translations of the structure, first the A-site sequence, then shifts the interstitials around to match
                        if accept == True:
                            comp_check=0
                            tmpB=temp[:,1].copy()
                            zero =0
                            one  =0
                            two  =0
                            three=0
                            for y in range(len(tmpA)):
                                if tmpC[y] in  a_2_dict[0]:
                                    zero = zero +1
                                if tmpC[y] in  a_2_dict[1]:
                                    one = one +1
                                if tmpC[y] in  a_2_dict[2]:
                                    two = two +1
                                if tmpC[y] in a_2_dict[3]:
                                    three=three+1
                            counter=list([zero,one,two])
                            if float(max(counter)) > float(len(tmpC)/2):
                                comp_check=1
                                
                            if float(max(counter)) <= float(len(tmpC)/2):                                
                                tmp_perm=tmpC
                                for i in range(len(tmp_perm)):### this bit checks the neighbours to see if they are of the same translation, would suggest that alter this so that if the are, then the translations are corrected                                    
                                    ### determine the translation of layer and neighbours
                                    loc1 = tmp_perm[i]
                                    try:
                                        loc2 = tmp_perm[i+1]
                                    except(IndexError):
                                        loc2 = tmp_perm[0]
                                    loc3=tmp_perm[i-1]
                                    idx1a=0
                                    idx2a=0
                                    idx3a=0
                                    neighbours=[]
                                    
                                    for j in range(len(list(a_2_dict.keys()))):
                                        if loc1 in a_2_dict[j]:
                                            idx1a=j
                                        if loc2 in a_2_dict[j]:
                                            idx2a=j
                                        if loc3 in a_2_dict[j]:
                                            idx3a=j
                                        if len(a_2_dict[j])!=0:
                                            neighbours.append(j)    
                                    if idx1a == idx2a: #if layer has equal translation to next layer in sequence, reject structure
                                        accept=False
                                        comp_check=1
                                        #try: ### was a routine to determine what the correct translation should be and correct for it, however, using this ended up making the code run slower
                                        #    neighbours.remove(idx2a)
                                        #    neighbours.remove(idx3a)
                                        #except(ValueError):
                                        #    print "3a",idx3a, "2a",idx2a
                                        #    sys.exit()
                                        #for j in range(len(a_dict.keys())):
                                        #    if loc1 in a_dict[j]:
                                        #        idx1b=j
                                        #tmpC[j]=a_dict[idx1b][idx1a]
                            tmpB=temp[:,1].copy()
                            ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
                            if comp_check !=1:
                                for b in range(len(tmpB)):
                                    trans=[0,0]
                                    ### identifying the A layer level with the B module
                                    if tmpC[b] in b_2_dict[0]:
                                        trans[0]=0
                                    if tmpC[b] in b_2_dict[1]:
                                        trans[0]=1
                                    if tmpC[b] in b_2_dict[2]:
                                        trans[0]=2
                                    if tmpC[b] in b_2_dict[3]:
                                        trans[0]=3
                                    ### identifying the neihbouring A module
                                    if b != max(range(len(tmpB))):
                                        if tmpC[b+1] in a_2_dict[0]:
                                            trans[1]=0
                                        if tmpC[b+1] in a_2_dict[1]:
                                            trans[1]=1
                                        if tmpC[b+1] in a_2_dict[2]:
                                            trans[1]=2
                                        if tmpC[b+1] in a_2_dict[3]:
                                            trans[1]=3
                                    if b+1 == len(tmpB):
                                        if tmpC[0] in a_2_dict[0]:
                                            trans[1]=0
                                        if tmpC[0] in a_2_dict[1]:
                                            trans[1]=1
                                        if tmpC[0] in a_2_dict[2]:
                                            trans[1]=2
                                        if tmpC[0] in a_2_dict[3]:
                                            trans[1]=3
                                    trans=numpy.array(trans)
                                    b_mod_idx=0
                                    if tmpB[b] in b_2_dict[0]:
                                        b_mod_idx=0
                                    if tmpB[b] in b_2_dict[1]:
                                        b_mod_idx=1
                                    if tmpB[b] in b_2_dict[2]:
                                        b_mod_idx=2
                                    if tmpB[b] in b_2_dict[3]:
                                        b_mod_idx=3
                                        
                                    ### sets the translations of the interstitial sites ###
                                    if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all():
                                        trans_1=2
                                        if trans_1 == b_mod_idx or b_mod_idx==3:
                                            pass
                                        if trans_1 != b_mod_idx and b_mod_idx != 3:
                                            for x in range(len(list(b_dict.keys()))):
                                                if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                                    loc1=x
                                                    loc2=b_dict[x]
                                                    tmpB[b]=b_dict[x][trans_1]
                                    if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                                        trans_1=1
                                        if trans_1 == b_mod_idx or b_mod_idx==3:
                                            pass
                                        if trans_1 != b_mod_idx and b_mod_idx != 3:
                                            for x in range(len(list(b_dict.keys()))):
                                                if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                                    loc1=x
                                                    loc2=b_dict[x]
                                                    tmpB[b]=b_dict[x][trans_1]
                                    if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                                        trans_1=0
                                        if trans_1 == b_mod_idx or b_mod_idx==3:
                                            pass
                                        if trans_1 != b_mod_idx and b_mod_idx != 3:
                                            for x in range(len(list(b_dict.keys()))):
                                                if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                                    loc1=x
                                                    loc2=b_dict[x]
                                                    tmpB[b]=b_dict[x][trans_1]
                                                    
                        if accept ==True:
                            for x in range(len(temp)):
                                temp[x,0]=tmpC[x]
                                temp[x,1]=tmpB[x]
                                struct=temp.copy()
                        #print struct
                        #sys.exit()
                        trial2=trial2+1
                    
                else:
                    pass
                
            
            else:        ### flag... in an ideal world... this is what we want all of the code to look like! (in terms of the spacings and commenting!)
                temp=struct.copy()
                SL = len(temp)
                if len(AMods) == 1:
                    ab = 1
                if len(BMods) == 1:
                    ab = 0
                else:
                    ab=choice(list(range(2)))
                tmp1=choice(list(range(len(temp[:,ab]))))
                tmp2=choice(list(range(len(temp[:,ab]))))
                                    
                l = 0
                while temp[tmp1,ab] == temp[tmp2,ab] and l <= 500:
                    tmp2=choice(list(range(len(temp[:,ab]))))
                    l = l+1

                    
                tmp1a=temp[tmp1,ab].copy()
                tmp2a=temp[tmp2,ab].copy()
                temp[tmp1,ab]=tmp2a
                temp[tmp2,ab]=tmp1a
                struct=temp.copy()
                
                
                
        if p_t == 'T2': ### type 2: the stacking order of the current structure should be randomised, retaining the current number and type of modules
            old_sl=SL
            temp = struct.copy()
            SL = len(temp)
            tmpA = temp[:,0].copy()
            tmpB = temp[:,1].copy()
            numpy.random.shuffle(tmpA)
            numpy.random.shuffle(tmpB)
            if types[t]=='hexagonal':
                zero =0
                one  =0
                two  =0
                three=0
                comp_check=0
                for y in range(len(tmpA)):
                    if tmpA[y] in  a_2_dict[0]:
                        zero = zero +1
                    if tmpA[y] in  a_2_dict[1]:
                        one = one +1
                    if tmpA[y] in  a_2_dict[2]:
                        two = two +1
                    if tmpA[y] in a_2_dict[3]:
                        three=three+1
                counter=list([zero,one,two])
                if float(max(counter)) > float(len(tmpA)/2):
                    comp_check=1
                if float(max(counter)) <= float(len(tmpA)/2):
                    b = 0
                    while b == 0:
                        b = 1
                        la=list(a_dict.values())
                        lb=list(b_dict.values())
                        numpy.random.shuffle(tmpA)
                        tmp_perm=tmpA
                        for i in range(len(tmp_perm)):
                            neighbours=[]
                            try:
                                loc1 = tmp_perm[i]
                                loc2 = tmp_perm[i+1]
                            except(IndexError):
                                loc1 = tmp_perm[i]
                                loc2 = tmp_perm[0]
                            neighbours.append(loc2)
                            loc3=''
                            for z in range(len(la)):
                                if loc2 in la[z]:
                                    loc3=la[z].index(loc2)
                                
                            for x in range(len(la)):
                                if len(la[x]) >> 0:
                                    if la[x][loc3] not in neighbours:
                                        neighbours.append(la[x][loc3])
                            
                            if loc1 in neighbours:
                                b = 0
                ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
                if comp_check !=1:
                    for b in range(len(tmpB)):
                        trans=[0,0]
                        ### identifying the A layer level with the B module
                        if tmpA[b] in b_2_dict[0]:
                            trans[0]=0
                        if tmpA[b] in b_2_dict[1]:
                            trans[0]=1
                        if tmpA[b] in b_2_dict[2]:
                            trans[0]=2
                        if tmpA[b] in b_2_dict[3]:
                            trans[0]=3
                        ### identifying the neihbouring a module
                        if b != max(range(len(tmpB))):
                            if tmpA[b+1] in a_2_dict[0]:
                                trans[1]=0
                            if tmpA[b+1] in a_2_dict[1]:
                                trans[1]=1
                            if tmpA[b+1] in a_2_dict[2]:
                                trans[1]=2
                            if tmpA[b+1] in a_2_dict[3]:
                                trans[1]=3
                        if b+1 == len(tmpB):
                            if tmpA[0] in a_2_dict[0]:
                                trans[1]=0
                            if tmpA[0] in a_2_dict[1]:
                                trans[1]=1
                            if tmpA[0] in a_2_dict[2]:
                                trans[1]=2
                            if tmpA[0] in a_2_dict[3]:
                                trans[1]=3
                        trans=numpy.array(trans)
                        b_mod_idx=0
                        if tmpB[b] in b_2_dict[0]:
                            b_mod_idx=0
                        if tmpB[b] in b_2_dict[1]:
                            b_mod_idx=1
                        if tmpB[b] in b_2_dict[2]:
                            b_mod_idx=2
                        if tmpB[b] in b_2_dict[3]:
                            b_mod_idx=3


                        if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all() == True:
                            trans_1=2
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                        if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                            trans_1=1
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                        if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                            trans_1=0
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                if comp_check !=1:
                    for x in range(len(struct)):
                        temp[x,0]=tmpA[x]
                        temp[x,1]=tmpB[x]
                    struct=temp.copy()
                
            
            else:
                numpy.random.shuffle(tmpA)
                numpy.random.shuffle(tmpB)
                temp[:,0]=tmpA.copy()
                temp[:,1]=tmpB.copy()
                struct=temp.copy()
            
        if p_t == 'T3': ### type 3: one layer in the current structure should be swapped for another, of a different type note that it should allow for the structure to maintain charge neutrality
            old_sl = SL
            Acomp_count = old_A_comp_count.copy()
            Bcomp_count = old_B_comp_count.copy()
            comp_check = 0
            if types[t]=='hexagonal':
                temp=struct.copy() ### creating copy of current structure
                ### determining how many non-vacency modules are in the run input
                lena=0
                lenb=0
                for i in range(len(list(a_dict.keys()))):
                    if len(a_dict[i]) >> 1:
                        lena=lena+1
                for i in range(len(list(b_dict.keys()))):
                    if len(b_dict[i]) >> 1:
                        lenb=lenb+1
                if lena == 1 and lenb == 1: ### when running with all of the other swap types on, this should then mean that (eventually) a different swap type is chosen if there is only one of each module type that is non-vacency
                    continue
                ab=[]
                if lena >>1 and lenb ==1: ### means new module should come from A-sites
                    ab=[0]
                if lenb >>1: ### means new module can come from with A-sites or interstitial sites, note that it has to have more than one species in the input file in order to swap, since the translations are determined by the A site layers
                    ab=[0,1]
                ab = choice(ab)
                swap_idx=choice(list(range(len(temp)))) ## chooses a layer to swap
                if ab ==0:### for swapping an A module
                    loc1a = temp[swap_idx-1,ab]#working out what the neighbouring modules are to the module chosen for swapping
                    try:
                        loc2a=temp[swap_idx+1,ab]
                    except(IndexError):
                        loc2a=temp[0,ab]
                    for i in range(len(list(a_2_dict.keys()))):### this determines the translations of the two neighbours
                        if loc1a in a_2_dict[i]:
                            loc1b=i
                        if loc2a in a_2_dict[i]:
                            loc2b=i
                    if loc1b != loc2b: ## if neighbouring layers are different translations, then the new module must be of the same translation
                        for i in range(len(list(a_2_dict.keys()))):
                            if temp[swap_idx,ab] in a_2_dict[i]:
                                loc3=i
                        swap_mod = choice(a_2_dict[loc3])
                        
                        trial = 0
                        while temp[swap_idx,ab]==swap_mod and trial <= 500:
                            swap_mod = choice(a_2_dict[loc3]) ### swaps the A module to a A module from the input file, inserting a new module of the same translation, since the translation is fixed by the neighbouring A-site layers
                            trial = trial +1
                        #    if trial == 500:
                        #        continue
                    
                    if loc1b == loc2b:
                        trial = list(range(3))
                        trial.remove(loc1b)
                        #trial.remove(loc2b)
                        loc3a=choice(trial)
                        swap_mod=choice(a_2_dict[loc3a])
                        trial2=0
                        while swap_mod==temp[swap_idx,ab] and trial2 <= 500:    
                            swap_mod=choice(a_2_dict[loc3a])
                            trial2 = trial2+1
                        #    if trial2 == 500:
                        #        continue
                    temp[swap_idx,ab]=swap_mod
                    
                    
                if ab ==1:
                    for i in range(len(list(b_2_dict.keys()))): ### determines the translation of the interstitial site layer
                        if temp[swap_idx,ab] in b_2_dict[i]:
                            loc1=i
                    swap_mod = choice(b_2_dict[loc1])
                    trial = 0
                    while temp[swap_idx,ab]==swap_mod and trial <= 500:
                        swap_mod = choice(b_2_dict[loc1]) ### swaps the interstitial module to a different interstitial from the input file, inserting a new module of the same translation, since the translation is fixed by the corrisponding A-site layers
                        trial = trial +1
                        #if trial == 500:
                        #    continue
                    temp[swap_idx,ab]=swap_mod
                    
            
            else: ### rountine for creating new structure for normal stacking types (i.e cubic)
                temp=struct.copy()
                SL = len(temp) 
                if len(AMods) == 1:
                    ab = 1
                    tmpA = list(range(nB))
                if len(BMods) == 1:
                    ab = 0
                    tmpA = list(range(nA))
                if len(BMods) and len(AMods) != 1:
                    ab=choice(list(range(2)))
                    if ab ==1:
                        tmpA=list(range(nB))
                    if ab ==0:
                        tmpA=list(range(nA))
                a=temp[:,ab]
                swap=choice(list(range(len(a))))
                swap_mod=temp[swap,ab]
                swap_new=swap_mod
                trial = 0
                while swap_new==swap_mod and trial <= 500:
                    swap_new=choice(tmpA)
                    trial = trial+1
                temp[swap,ab]=swap_new
            ### this is the composition checking part of the routine
            keys=[]
            gen_comp=[]
            for i in range(len(temp)):
                tmp=Acomp_count[stack][temp[i,0]]
                for x in range(len(Acomp_count[stack][temp[i,0]])):
                    gen_comp.append(Acomp_count[stack][temp[i,0]][x])
                for x in range(len(tmp)):
                    if tmp[x] not in keys:
                        keys.append(tmp[x])
                #print stack
                #print temp[i,1]
                #print Bcomp_count
                tmp=Bcomp_count[stack][temp[i,1]]
                
                for x in range(len(Bcomp_count[stack][temp[i,1]])):
                    gen_comp.append(Bcomp_count[stack][temp[i,1]][x])
                for x in range(len(tmp)):
                    if tmp[x] not in keys:
                        keys.append(tmp[x])
            trial_comp={}
            for i in range(len(keys)):
                tmpcount=gen_comp.count(keys[i])
                trial_comp[keys[i]]=tmpcount
            min_symb=min(trial_comp,key=trial_comp.get)
            norm = trial_comp[min_symb]
            for sym in range(len(keys)):
                tmp=float(float(trial_comp[keys[sym]])/float(norm))
                trial_comp[keys[sym]]=tmp
            b = trial_comp==composition    
            if b == True: 
                if types[t] == 'hexagonal':
                    if comp_check != 1:
                        struct=temp.copy()
                else:
                    struct=temp.copy()
                
        if p_t == 'T4': ### type 4: generate an entire new set of layers, maintaining the current stacking length
            a = False
            SL = len(struct)
            old_sl = SL        
            z = 1
            accept =False
            make = False
            tmax = tmax
            #old_A_comp_count = Acomp_count.copy()
            #old_B_comp_count = Bcomp_count.copy()
            Acomp_count=old_A_comp_count[types[t]]
            Bcomp_count=old_B_comp_count[types[t]]
            n = 0
            while make==False:
                z = 0
                while z < tmax:
                    struct=make_random_structure(AMods,BMods,composition,[SL])
                    if type(struct) != numpy.ndarray:
                        z+=1
                    if type(struct) == numpy.ndarray:
                        a=True
                        make = True
                        if types[t] != 'hexagonal':
                            accept=True

                        break
                if z == tmax:
                    continue                    
            
                #generate first structure at random, accecpting the first generated structure which is charge neutral#
                #when the structure type is hexagonal, now need to go through and setup the second array defining the translations
                if types[t] == 'hexagonal':
                    comp_check=0
                    struct_Amods=[]
                    struct_Bmods=[]
                    for mod in range(len(struct)):
                        struct_Amods.append(struct[mod][0])
                        struct_Bmods.append(struct[mod][1])
                    mods = []
                    ModA=()
                    ModB=()
                    for i in range(len(A_Mods[types[t]])):
                            ModA += (struct_Amods.count(i),)
                    for i in range(len(B_Mods[types[t]])):
                            ModB += (struct_Bmods.count(i),)
                    mods=[ModA,ModB]
                    tmpA=struct_Amods
                    tmpB=struct_Bmods
                    Aset=mods[0]
                    Bset=mods[1]
                    zero =0
                    one  =0
                    two  =0
                    three=0
                    for y in range(len(tmpA)):
                        if tmpA[y] in  a_2_dict[0]:
                            zero = zero +1
                        if tmpA[y] in  a_2_dict[1]:
                            one = one +1
                        if tmpA[y] in  a_2_dict[2]:
                            two = two +1
                        if tmpA[y] in a_2_dict[3]:
                            three=three+1
                    test=list([zero,one,two])
                    if float(max(test)) > float(len(tmpA)/2):
                        comp_check=1
                    if float(max(test)) <= float(len(tmpA)/2):
                        b = 0
                        while b == 0:
                            b = 1
                            numpy.random.shuffle(tmpA)
                            la=list(a_dict.values())
                            lb=list(b_dict.values())
                            tmp_perm=tmpA
                            for i in range(len(tmp_perm)):
                                neighbours=[]
                                try:
                                    loc1 = tmp_perm[i]
                                    loc2 = tmp_perm[i+1]
                                except(IndexError):
                                    loc1 = tmp_perm[i]
                                    loc2 = tmp_perm[0]
                                neighbours.append(loc2)
                                loc3=''
                                for z in range(len(la)):
                                    if loc2 in la[z]:
                                        loc3=la[z].index(loc2)
                                    
                                for x in range(len(la)):
                                    if len(la[x]) >> 0:
                                        if la[x][loc3] not in neighbours:
                                            neighbours.append(la[x][loc3])
                                   
                                if loc1 in neighbours:
                                    b = 0
                    ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
                    if comp_check !=1:
                        numpy.random.shuffle(tmpB)
                        for b in range(len(tmpB)):
                            trans=[0,0]
                            ### identifying the A layer level with the B module
                            if tmpA[b] in b_2_dict[0]:
                                trans[0]=0
                            if tmpA[b] in b_2_dict[1]:
                                trans[0]=1
                            if tmpA[b] in b_2_dict[2]:
                                trans[0]=2
                            if tmpA[b] in b_2_dict[3]:
                                trans[0]=3
                            ### identifying the neihbouring a module
                            if b != max(range(len(tmpB))):
                                if tmpA[b+1] in a_2_dict[0]:
                                    trans[1]=0
                                if tmpA[b+1] in a_2_dict[1]:
                                    trans[1]=1
                                if tmpA[b+1] in a_2_dict[2]:
                                    trans[1]=2
                                if tmpA[b+1] in a_2_dict[3]:
                                    trans[1]=3
                            if b+1 == len(tmpB):
                                if tmpA[0] in a_2_dict[0]:
                                    trans[1]=0
                                if tmpA[0] in a_2_dict[1]:
                                    trans[1]=1
                                if tmpA[0] in a_2_dict[2]:
                                    trans[1]=2
                                if tmpA[0] in a_2_dict[3]:
                                    trans[1]=3
                            trans=numpy.array(trans)
                            b_mod_idx=0
                            if tmpB[b] in b_2_dict[0]:
                                b_mod_idx=0
                            if tmpB[b] in b_2_dict[1]:
                                b_mod_idx=1
                            if tmpB[b] in b_2_dict[2]:
                                b_mod_idx=2
                            if tmpB[b] in b_2_dict[3]:
                                b_mod_idx=3
                
                
                            if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all() == True:
                                trans_1=2
                                if trans_1 == b_mod_idx or b_mod_idx==3:
                                    pass
                                if trans_1 != b_mod_idx and b_mod_idx != 3:
                                    for x in range(len(list(b_dict.keys()))):
                                        if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                            loc1=x
                                            loc2=b_dict[x]
                                            tmpB[b]=b_dict[x][trans_1]
                            if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                                trans_1=1
                                if trans_1 == b_mod_idx or b_mod_idx==3:
                                    pass
                                if trans_1 != b_mod_idx and b_mod_idx != 3:
                                    for x in range(len(list(b_dict.keys()))):
                                        if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                            loc1=x
                                            loc2=b_dict[x]
                                            tmpB[b]=b_dict[x][trans_1]
                            if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                                trans_1=0
                                if trans_1 == b_mod_idx or b_mod_idx==3:
                                    pass
                                if trans_1 != b_mod_idx and b_mod_idx != 3:
                                    for x in range(len(list(b_dict.keys()))):
                                        if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                            loc1=x
                                            loc2=b_dict[x]
                                            tmpB[b]=b_dict[x][trans_1]
                    if comp_check !=1:
                        for x in range(len(struct)):
                            struct[x,0]=tmpA[x]
                            struct[x,1]=tmpB[x]
                        ### ok, now correctly identifying what the b-site translations should be and what they are.... now need to go back and find a way if they are different to change the module for a diffreent translation, but of the same chemical composition... think this should be done using the b_dict dictionary, so you can pull a value out of the same key, but this time with the correct index
                
                    else:
                        tmpA=[]
                        tmpB=[]
                        struct_modA=()
                        struct_modB=()
                        for i in range(len(struct)):
                            tmpA.append(struct[i][0])
                            tmpB.append(struct[i][1])
                        for i in range(len(A_Mods)):
                            struct_modA+=(tmpA.count(i),)
                        for i in range(len(B_Mods)):
                            struct_modB+=(tmpB.count(i),)
                        mods=(tmpA,tmpB)
                        Aset=mods[0]
                        Bset=mods[1]
                        tmpA=numpy.array(tmpA)
                        tmpB=numpy.array(tmpB)
                        struct_charge=[0] #### is this due to me removing the charge counting? as I can't seem to find it now... ####
                        struct_charge=numpy.array(struct_charge)
                        comp_check=0
                    
                    ### this is the composition checking section, so this should apply to all structure types ###
                    keys=[]
                    gen_comp=[]                
                    Acomp_count=old_A_comp_count.copy()
                    Bcomp_count=old_B_comp_count.copy()
                    #sys.exit()
                    for i in range(len(struct)):
                        tmp=Acomp_count[stack][struct[i,0]]
                        for x in range(len(Acomp_count[stack][struct[i,0]])):
                            gen_comp.append(Acomp_count[stack][struct[i,0]][x])
                        for x in range(len(tmp)):
                            if tmp[x] not in keys:
                                keys.append(tmp[x])
                        tmp=Bcomp_count[stack][struct[i,1]]
                        for x in range(len(Bcomp_count[stack][struct[i,1]])):
                            gen_comp.append(Bcomp_count[stack][struct[i,1]][x])
                        for x in range(len(tmp)):
                            if tmp[x] not in keys:
                                keys.append(tmp[x])
                    trial_comp={}
                    for i in range(len(keys)):
                        tmpcount=gen_comp.count(keys[i])
                        trial_comp[keys[i]]=tmpcount
                    min_symb=min(trial_comp,key=trial_comp.get)
                    norm = trial_comp[min_symb]
                    for sym in range(len(keys)):
                        tmp=float(float(trial_comp[keys[sym]])/float(norm))
                        trial_comp[keys[sym]]=tmp
                        
                    a=trial_comp==composition
                    if a==False:
                        comp_check=1
                    
                    if comp_check==0:
                        test =  0
                        if str(s_type) == 'hexagonal':
                            test = check_struct_hexagonal(struct,perms,accepted=[],a_dict=a_dict,b_dict=b_dict,A_MODS=A_Mods['hexagonal'],B_MODS=B_Mods['hexagonal'],ABC=ABC)
                        if test == 0:
                            new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
                            if min_bond_num == '':
                                distances=get_distances(new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist))
                                if min(distances) <= min_bond:
                                    test = 1
                            if min_bond_num != '':
### min_bond_nu,     should be a parameter to use to specifiy a minimum number of bonds for different species, expressed as a dictionary. e.g. min_bond_num = {"Al":[4,'O',3]}, note that this specifies min of 4 bonds for Al to Oxygen, with a cufoff of 3 Ang.
                                if len(list(min_bond_num.keys())) >> 0:
                                    bonds=[]
                                    for i in range(len(list(min_bond_num.keys()))):
                                        cutoff=float(min_bond_num[list(min_bond_num.keys())[i]][2])
                                        num=min_bond_num[list(min_bond_num.keys())[i]][0]
                                        s1=list(min_bond_num.keys())[i]
                                        s2=min_bond_num[list(min_bond_num.keys())[i]][1]
                                        s1a=[]
                                        s2a=[]
                                        for a in range(len(new_atoms)):
                                            if new_atoms[a].symbol==s1:
                                                s1a.append(a)
                                            if new_atoms[a].symbol==s2:
                                                s2a.append(a)
                                        for x in range(len(s1a)):
                                            rij=[]
                                            for y in range(len(s2a)):
                                                tmp=new_atoms.get_distance(s1a[x],s2a[y],mic=True)
                                                if tmp <= cutoff and tmp >= min_bond:
                                                    rij.append(tmp)
                                            bonds.append(rij)
                                    for x in range(len(bonds)):
                                        #print len(bonds[x])
                                        if len(bonds[x]) < num:
                                            test=1
                                    
                        if test == 0:
                            z=0
                            accept = True
                            make = True
                        else:
                            z=1
                            make = False
                n = n+1
                #print n
                if n == 500000:
                    break
        if p_t == 'T5': ### type 5: multiply the structure along c, retaining the original structure, 50:50 chance of tiling or adding in new modules
            a = False
            SL = len(struct)
            old_sl=len(struct)
            if len(sl) == 1:
                force = delay
                f_type='low'
                continue
            # work out if it's possible to elongate the structure within the restrictions of the system
            possible_multipliers=[]
            for i in sl:
            	if i % SL == 0:
            		possible_multipliers.append(i)
            possible_multipliers.remove(SL)
            if len(possible_multipliers) == 0: # if cannot extend strcture, default to new random structure
            	p_t = 'T6'
            	continue
        
            new_length=choice(possible_multipliers)
            
            setting=choice([1,2]) # 1 tile current structure, 2 append current structure with random new modules
            setting = 2
            mult=int(new_length/SL)
            
            if setting == 2:
                dummy_sl=[int((mult-1)*old_sl)]
                #print("dummy sl: \n",dummy_sl)
                make = False
                first_section=struct.copy()
                #print(AMods)
                #print(BMods)
                #print(composition)
                #print(dummy_sl)
                while make==False:
                    z = 0
                    while z < tmax:
                        #print(z)
                        second_section=make_random_structure(AMods,BMods,composition,dummy_sl)
                        #print(second_section)
                        #print(type(second_section)==numpy.ndarray)
                        
                        if type(second_section) != numpy.ndarray:
                            z+=1
                        if type(second_section) == numpy.ndarray:
                            #print("hello")
                            make=True
                            break
                        
                        #sys.exit()
                        
                    if z == tmax:
                    	  setting = 1
                    	  break
                #print("first section: \n",first_section)
                #print("second_section: \n",second_section)
                new_struct=list(first_section.copy())
                for i in second_section:
                	 new_struct.append(i)
                new_struct=numpy.array(new_struct)
                #print(new_struct)
                #sys.exit()
                
            if setting == 1:
            	#print("starting structure: \n",struct)
            	#print("mutliplier: \n",mult)
            	mods_to_tile=struct.copy()
            	new_struct=list(struct.copy())
            	for i in range(mult-1):
            		for j in mods_to_tile:
            			new_struct.append(j)
            	new_struct=numpy.array(new_struct)
            	#print("new_structure: \n",new_struct)
            	#print("new length :\n",len(new_struct))
            
            
            if types[t] != 'hexagonal':
            	a=True
            	make=True
            	accept=True
            	SL = len(new_struct)
            	struct=new_struct.copy()
            	
            if types[t] == 'hexagonal':
                comp_check=0
                struct_Amods=[]
                struct_Bmods=[]
                for mod in range(len(struct)):
                    struct_Amods.append(struct[mod][0])
                    struct_Bmods.append(struct[mod][1])
                mods = []
                ModA=()
                ModB=()
                for i in range(len(A_Mods[types[t]])):
                        ModA += (struct_Amods.count(i),)
                for i in range(len(B_Mods[types[t]])):
                        ModB += (struct_Bmods.count(i),)
                mods=[ModA,ModB]
                tmpA=struct_Amods
                tmpB=struct_Bmods
                Aset=mods[0]
                Bset=mods[1]
                zero =0
                one  =0
                two  =0
                three=0
                for y in range(len(tmpA)):
                    if tmpA[y] in  a_2_dict[0]:
                        zero = zero +1
                    if tmpA[y] in  a_2_dict[1]:
                        one = one +1
                    if tmpA[y] in  a_2_dict[2]:
                        two = two +1
                    if tmpA[y] in a_2_dict[3]:
                        three=three+1
                test=list([zero,one,two])
                if float(max(test)) > float(len(tmpA)/2):
                    comp_check=1
                if float(max(test)) <= float(len(tmpA)/2):
                    b = 0
                    while b == 0:
                        b = 1
                        numpy.random.shuffle(tmpA)
                        la=list(a_dict.values())
                        lb=list(b_dict.values())
                        tmp_perm=tmpA
                        for i in range(len(tmp_perm)):
                            neighbours=[]
                            try:
                                loc1 = tmp_perm[i]
                                loc2 = tmp_perm[i+1]
                            except(IndexError):
                                loc1 = tmp_perm[i]
                                loc2 = tmp_perm[0]
                            neighbours.append(loc2)
                            loc3=''
                            for z in range(len(la)):
                                if loc2 in la[z]:
                                    loc3=la[z].index(loc2)
                                
                            for x in range(len(la)):
                                if len(la[x]) >> 0:
                                    if la[x][loc3] not in neighbours:
                                        neighbours.append(la[x][loc3])
                               
                            if loc1 in neighbours:
                                b = 0
                ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
                if comp_check !=1:
                    numpy.random.shuffle(tmpB)
                    for b in range(len(tmpB)):
                        trans=[0,0]
                        ### identifying the A layer level with the B module
                        if tmpA[b] in b_2_dict[0]:
                            trans[0]=0
                        if tmpA[b] in b_2_dict[1]:
                            trans[0]=1
                        if tmpA[b] in b_2_dict[2]:
                            trans[0]=2
                        if tmpA[b] in b_2_dict[3]:
                            trans[0]=3
                        ### identifying the neihbouring a module
                        if b != max(range(len(tmpB))):
                            if tmpA[b+1] in a_2_dict[0]:
                                trans[1]=0
                            if tmpA[b+1] in a_2_dict[1]:
                                trans[1]=1
                            if tmpA[b+1] in a_2_dict[2]:
                                trans[1]=2
                            if tmpA[b+1] in a_2_dict[3]:
                                trans[1]=3
                        if b+1 == len(tmpB):
                            if tmpA[0] in a_2_dict[0]:
                                trans[1]=0
                            if tmpA[0] in a_2_dict[1]:
                                trans[1]=1
                            if tmpA[0] in a_2_dict[2]:
                                trans[1]=2
                            if tmpA[0] in a_2_dict[3]:
                                trans[1]=3
                        trans=numpy.array(trans)
                        b_mod_idx=0
                        if tmpB[b] in b_2_dict[0]:
                            b_mod_idx=0
                        if tmpB[b] in b_2_dict[1]:
                            b_mod_idx=1
                        if tmpB[b] in b_2_dict[2]:
                            b_mod_idx=2
                        if tmpB[b] in b_2_dict[3]:
                            b_mod_idx=3
            
            
                        if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all() == True:
                            trans_1=2
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                        if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                            trans_1=1
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                        if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                            trans_1=0
                            if trans_1 == b_mod_idx or b_mod_idx==3:
                                pass
                            if trans_1 != b_mod_idx and b_mod_idx != 3:
                                for x in range(len(list(b_dict.keys()))):
                                    if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                        loc1=x
                                        loc2=b_dict[x]
                                        tmpB[b]=b_dict[x][trans_1]
                if comp_check !=1:
                    for x in range(len(struct)):
                        struct[x,0]=tmpA[x]
                        struct[x,1]=tmpB[x]
                    ### ok, now correctly identifying what the b-site translations should be and what they are.... now need to go back and find a way if they are different to change the module for a diffreent translation, but of the same chemical composition... think this should be done using the b_dict dictionary, so you can pull a value out of the same key, but this time with the correct index
            
                else:
                    tmpA=[]
                    tmpB=[]
                    struct_modA=()
                    struct_modB=()
                    for i in range(len(struct)):
                        tmpA.append(struct[i][0])
                        tmpB.append(struct[i][1])
                    for i in range(len(A_Mods)):
                        struct_modA+=(tmpA.count(i),)
                    for i in range(len(B_Mods)):
                        struct_modB+=(tmpB.count(i),)
                    mods=(tmpA,tmpB)
                    Aset=mods[0]
                    Bset=mods[1]
                    tmpA=numpy.array(tmpA)
                    tmpB=numpy.array(tmpB)
                    struct_charge=[0] #### is this due to me removing the charge counting? as I can't seem to find it now... ####
                    struct_charge=numpy.array(struct_charge)
                    comp_check=0
                
                ### this is the composition checking section, so this should apply to all structure types ###
                keys=[]
                gen_comp=[]                
                Acomp_count=old_A_comp_count.copy()
                Bcomp_count=old_B_comp_count.copy()
                #sys.exit()
                for i in range(len(struct)):
                    tmp=Acomp_count[stack][struct[i,0]]
                    for x in range(len(Acomp_count[stack][struct[i,0]])):
                        gen_comp.append(Acomp_count[stack][struct[i,0]][x])
                    for x in range(len(tmp)):
                        if tmp[x] not in keys:
                            keys.append(tmp[x])
                    tmp=Bcomp_count[stack][struct[i,1]]
                    for x in range(len(Bcomp_count[stack][struct[i,1]])):
                        gen_comp.append(Bcomp_count[stack][struct[i,1]][x])
                    for x in range(len(tmp)):
                        if tmp[x] not in keys:
                            keys.append(tmp[x])
                trial_comp={}
                for i in range(len(keys)):
                    tmpcount=gen_comp.count(keys[i])
                    trial_comp[keys[i]]=tmpcount
                min_symb=min(trial_comp,key=trial_comp.get)
                norm = trial_comp[min_symb]
                for sym in range(len(keys)):
                    tmp=float(float(trial_comp[keys[sym]])/float(norm))
                    trial_comp[keys[sym]]=tmp
                    
                a=trial_comp==composition
                if a==False:
                    comp_check=1
                
                if comp_check==0:
                    test =  0
                    if str(s_type) == 'hexagonal':
                        test = check_struct_hexagonal(struct,perms,accepted=[],a_dict=a_dict,b_dict=b_dict,A_MODS=A_Mods['hexagonal'],B_MODS=B_Mods['hexagonal'],ABC=ABC)
                    if test == 0:
                        new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
                        if min_bond_num == '':
                            distances=get_distances(new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist))
                            if min(distances) <= min_bond:
                                test = 1
                        if min_bond_num != '':
### min_bond     should be a parameter to use to specifiy a minimum number of bonds for different species, expressed as a dictionary. e.g. min_bond_num = {"Al":[4,'O',3]}, note that this specifies min of 4 bonds for Al to Oxygen, with a cufoff of 3 Ang.
                            if len(list(min_bond_num.keys())) >> 0:
                                bonds=[]
                                for i in range(len(list(min_bond_num.keys()))):
                                    cutoff=float(min_bond_num[list(min_bond_num.keys())[i]][2])
                                    num=min_bond_num[list(min_bond_num.keys())[i]][0]
                                    s1=list(min_bond_num.keys())[i]
                                    s2=min_bond_num[list(min_bond_num.keys())[i]][1]
                                    s1a=[]
                                    s2a=[]
                                    for a in range(len(new_atoms)):
                                        if new_atoms[a].symbol==s1:
                                            s1a.append(a)
                                        if new_atoms[a].symbol==s2:
                                            s2a.append(a)
                                    for x in range(len(s1a)):
                                        rij=[]
                                        for y in range(len(s2a)):
                                            tmp=new_atoms.get_distance(s1a[x],s2a[y],mic=True)
                                            if tmp <= cutoff and tmp >= min_bond:
                                                rij.append(tmp)
                                        bonds.append(rij)
                                for x in range(len(bonds)):
                                    #print len(bonds[x])
                                    if len(bonds[x]) < num:
                                        test=1
                                
                    if test == 0:
                        z=0
                        accept = True
                        make = True
                    else:
                        z=1
                        make = False            		
        
        
        if p_t == 'T6': ### type 6: change the stacking length, generating a new module set
            old_struct=struct[:]
            old_sl=len(old_struct)
            a = False
            SL = len(struct)
            if len(sl) == 1:
                force = delay
                f_type='low'
                continue
            old_sl = SL
            tmp_sl = sl.copy()
            tmp_sl.remove(old_sl)
            n_sl = choice(tmp_sl)
            accept =False
            make = False
            tmax = tmax
            Acomp_count=old_A_comp_count[types[t]]
            Bcomp_count=old_B_comp_count[types[t]]
            n=0
            #first section, choose a valid module set for the new structure
            while make==False:
                z = 0
                while z < tmax:
                    struct=make_random_structure(AMods,BMods,composition,tmp_sl)
                    if type(struct) != numpy.ndarray:
                        z+=1
                    if type(struct) == numpy.ndarray:
                        make=True
                        a=True
                        if types[t] != 'hexagonal':
                            accept=True
                        break
                if z == tmax:
                    continue  
				#sys.exit()

                #generate first structure at random, accecpting the first generated structure which is charge neutral#
                #when the structure type is hexagonal, now need to go through and setup the second array defining the translations
                if types[t] == 'hexagonal':
                    comp_check=0
                    struct_Amods=[]
                    struct_Bmods=[]
                    for mod in range(len(struct)):
                        struct_Amods.append(struct[mod][0])
                        struct_Bmods.append(struct[mod][1])
                    mods = []
                    ModA=()
                    ModB=()
                    for i in range(len(A_Mods[types[t]])):
                            ModA += (struct_Amods.count(i),)
                    for i in range(len(B_Mods[types[t]])):
                            ModB += (struct_Bmods.count(i),)
                    mods=[ModA,ModB]
                    tmpA=struct_Amods
                    tmpB=struct_Bmods
                    Aset=mods[0]
                    Bset=mods[1]
                    zero =0
                    one  =0
                    two  =0
                    three=0
                    for y in range(len(tmpA)):
                        if tmpA[y] in  a_2_dict[0]:
                            zero = zero +1
                        if tmpA[y] in  a_2_dict[1]:
                            one = one +1
                        if tmpA[y] in  a_2_dict[2]:
                            two = two +1
                        if tmpA[y] in a_2_dict[3]:
                            three=three+1
                    test=list([zero,one,two])
                    if float(max(test)) > float(len(tmpA)/2):
                        comp_check=1
                    if float(max(test)) <= float(len(tmpA)/2):
                        b = 0
                        while b == 0:
                            b = 1
                            numpy.random.shuffle(tmpA)
                            la=list(a_dict.values())
                            lb=list(b_dict.values())
                            tmp_perm=tmpA
                            for i in range(len(tmp_perm)):
                                neighbours=[]
                                try:
                                    loc1 = tmp_perm[i]
                                    loc2 = tmp_perm[i+1]
                                except(IndexError):
                                    loc1 = tmp_perm[i]
                                    loc2 = tmp_perm[0]
                                neighbours.append(loc2)
                                loc3=''
                                for z in range(len(la)):
                                    if loc2 in la[z]:
                                        loc3=la[z].index(loc2)
                                    
                                for x in range(len(la)):
                                    if len(la[x]) >> 0:
                                        if la[x][loc3] not in neighbours:
                                            neighbours.append(la[x][loc3])
                                   
                                if loc1 in neighbours:
                                    b = 0
                    ### now going through the generated B-metal sequence and checking that it is of the correct translation for the A-site sequence
                    if comp_check !=1:
                        numpy.random.shuffle(tmpB)
                        for b in range(len(tmpB)):
                            trans=[0,0]
                            ### identifying the A layer level with the B module
                            if tmpA[b] in b_2_dict[0]:
                                trans[0]=0
                            if tmpA[b] in b_2_dict[1]:
                                trans[0]=1
                            if tmpA[b] in b_2_dict[2]:
                                trans[0]=2
                            if tmpA[b] in b_2_dict[3]:
                                trans[0]=3
                            ### identifying the neihbouring a module
                            if b != max(range(len(tmpB))):
                                if tmpA[b+1] in a_2_dict[0]:
                                    trans[1]=0
                                if tmpA[b+1] in a_2_dict[1]:
                                    trans[1]=1
                                if tmpA[b+1] in a_2_dict[2]:
                                    trans[1]=2
                                if tmpA[b+1] in a_2_dict[3]:
                                    trans[1]=3
                            if b+1 == len(tmpB):
                                if tmpA[0] in a_2_dict[0]:
                                    trans[1]=0
                                if tmpA[0] in a_2_dict[1]:
                                    trans[1]=1
                                if tmpA[0] in a_2_dict[2]:
                                    trans[1]=2
                                if tmpA[0] in a_2_dict[3]:
                                    trans[1]=3
                            trans=numpy.array(trans)
                            b_mod_idx=0
                            if tmpB[b] in b_2_dict[0]:
                                b_mod_idx=0
                            if tmpB[b] in b_2_dict[1]:
                                b_mod_idx=1
                            if tmpB[b] in b_2_dict[2]:
                                b_mod_idx=2
                            if tmpB[b] in b_2_dict[3]:
                                b_mod_idx=3
                
                
                            if numpy.equal(trans,(numpy.array([0,1]))).all() or numpy.equal(trans,(numpy.array([1,0]))).all() == True:
                                trans_1=2
                                if trans_1 == b_mod_idx or b_mod_idx==3:
                                    pass
                                if trans_1 != b_mod_idx and b_mod_idx != 3:
                                    for x in range(len(list(b_dict.keys()))):
                                        if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                            loc1=x
                                            loc2=b_dict[x]
                                            tmpB[b]=b_dict[x][trans_1]
                            if numpy.equal(trans,(numpy.array([0,2]))).all() or numpy.equal(trans,(numpy.array([2,0]))).all() == True:
                                trans_1=1
                                if trans_1 == b_mod_idx or b_mod_idx==3:
                                    pass
                                if trans_1 != b_mod_idx and b_mod_idx != 3:
                                    for x in range(len(list(b_dict.keys()))):
                                        if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                            loc1=x
                                            loc2=b_dict[x]
                                            tmpB[b]=b_dict[x][trans_1]
                            if numpy.equal(trans,(numpy.array([1,2]))).all() or numpy.equal(trans,(numpy.array([2,1]))).all() == True:
                                trans_1=0
                                if trans_1 == b_mod_idx or b_mod_idx==3:
                                    pass
                                if trans_1 != b_mod_idx and b_mod_idx != 3:
                                    for x in range(len(list(b_dict.keys()))):
                                        if tmpB[b] in b_dict[list(b_dict.keys())[x]]:
                                            loc1=x
                                            loc2=b_dict[x]
                                            tmpB[b]=b_dict[x][trans_1]
                    if comp_check !=1:
                        for x in range(len(struct)):
                            struct[x,0]=tmpA[x]
                            struct[x,1]=tmpB[x]
                        ### ok, now correctly identifying what the b-site translations should be and what they are.... now need to go back and find a way if they are different to change the module for a diffreent translation, but of the same chemical composition... think this should be done using the b_dict dictionary, so you can pull a value out of the same key, but this time with the correct index
                
                    else:
                        tmpA=[]
                        tmpB=[]
                        struct_modA=()
                        struct_modB=()
                        for i in range(len(struct)):
                            tmpA.append(struct[i][0])
                            tmpB.append(struct[i][1])
                        for i in range(len(A_Mods)):
                            struct_modA+=(tmpA.count(i),)
                        for i in range(len(B_Mods)):
                            struct_modB+=(tmpB.count(i),)
                        mods=(tmpA,tmpB)
                        Aset=mods[0]
                        Bset=mods[1]
                        tmpA=numpy.array(tmpA)
                        tmpB=numpy.array(tmpB)
                        struct_charge=[0] #### is this due to me removing the charge counting? as I can't seem to find it now... ####
                        struct_charge=numpy.array(struct_charge)
                        comp_check=0
                    
                    ### this is the composition checking section, so this should apply to all structure types ###
                    keys=[]
                    gen_comp=[]                
                    Acomp_count=old_A_comp_count.copy()
                    Bcomp_count=old_B_comp_count.copy()
                    #sys.exit()
                    for i in range(len(struct)):
                        tmp=Acomp_count[stack][struct[i,0]]
                        for x in range(len(Acomp_count[stack][struct[i,0]])):
                            gen_comp.append(Acomp_count[stack][struct[i,0]][x])
                        for x in range(len(tmp)):
                            if tmp[x] not in keys:
                                keys.append(tmp[x])
                        tmp=Bcomp_count[stack][struct[i,1]]
                        for x in range(len(Bcomp_count[stack][struct[i,1]])):
                            gen_comp.append(Bcomp_count[stack][struct[i,1]][x])
                        for x in range(len(tmp)):
                            if tmp[x] not in keys:
                                keys.append(tmp[x])
                    trial_comp={}
                    for i in range(len(keys)):
                        tmpcount=gen_comp.count(keys[i])
                        trial_comp[keys[i]]=tmpcount
                    min_symb=min(trial_comp,key=trial_comp.get)
                    norm = trial_comp[min_symb]
                    for sym in range(len(keys)):
                        tmp=float(float(trial_comp[keys[sym]])/float(norm))
                        trial_comp[keys[sym]]=tmp
                        
                    a=trial_comp==composition
                    if a==False:
                        comp_check=1
                    
                    if comp_check==0:
                        test =  0
                        if str(s_type) == 'hexagonal':
                            test = check_struct_hexagonal(struct,perms,accepted=[],a_dict=a_dict,b_dict=b_dict,A_MODS=A_Mods['hexagonal'],B_MODS=B_Mods['hexagonal'],ABC=ABC)
                        if test == 0:
                            new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
                            if min_bond_num == '':
                                distances=get_distances(new_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist))
                                if min(distances) <= min_bond:
                                    test = 1
                            if min_bond_num != '':
### min_bond_nu,     should be a parameter to use to specifiy a minimum number of bonds for different species, expressed as a dictionary. e.g. min_bond_num = {"Al":[4,'O',3]}, note that this specifies min of 4 bonds for Al to Oxygen, with a cufoff of 3 Ang.
                                if len(list(min_bond_num.keys())) >> 0:
                                    bonds=[]
                                    for i in range(len(list(min_bond_num.keys()))):
                                        cutoff=float(min_bond_num[list(min_bond_num.keys())[i]][2])
                                        num=min_bond_num[list(min_bond_num.keys())[i]][0]
                                        s1=list(min_bond_num.keys())[i]
                                        s2=min_bond_num[list(min_bond_num.keys())[i]][1]
                                        s1a=[]
                                        s2a=[]
                                        for a in range(len(new_atoms)):
                                            if new_atoms[a].symbol==s1:
                                                s1a.append(a)
                                            if new_atoms[a].symbol==s2:
                                                s2a.append(a)
                                        for x in range(len(s1a)):
                                            rij=[]
                                            for y in range(len(s2a)):
                                                tmp=new_atoms.get_distance(s1a[x],s2a[y],mic=True)
                                                if tmp <= cutoff and tmp >= min_bond:
                                                    rij.append(tmp)
                                            bonds.append(rij)
                                    for x in range(len(bonds)):
                                        #print len(bonds[x])
                                        if len(bonds[x]) < num:
                                            test=1
                                    
                        if test == 0:
                            z=0
                            accept = True
                            make = True
                        else:
                            z=1
                            make = False
                n = n+1
                #print n
                if n == 500000:
                    break
                     

        p_log.append(p_t)
        ### note that there is at least another pertubation type possible here: type 7: to change the stacking type of the structure: i.e. to go from building a cubic based perovskite to a hexagonal one ###
        
        ###routine to check if new atoms has been generated before, when have implimented other permutation options, will impliment it such that if k >= 1000 p_t = 3,4 or 6 from above in the next cycle, and then only breaking the loop if this results in no new permutations, note that after types 5-7 have been implimented, will need to alter the way that check_cubic works, probably by storing previous permutations by stacking length, dictionary of lists, with the key(s) being the stacking lengths previously trialled?
        test = 1
        #print k
        if check =='T':
            if types[t] == 'hexagonal':
                #print "hexagonal"
                test = check_struct_hexagonal(struct,perms,accepted,a_dict=a_dict,b_dict=b_dict,A_MODS=A_Mods['hexagonal'],B_MODS=B_Mods['hexagonal'],ABC=ABC)
            else:
                test=check_struct_cubic(struct,perms,accepted)    
            if test == 0:
                #print "struct", struct
                SL=len(struct)
                trial_atoms= perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
                distances=get_distances(new_atoms=trial_atoms)
                if min_bond_num == '':
                    if min(distances) <= min_bond:
                        test = 1
                #if not min(distances) <= min_bond:    
                if min_bond_num != '':
                    trial_atoms = perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
                    #print "flag3"
                    if len(list(min_bond_num.keys())) >> 0:
                        bonds=[]
                        for i in range(len(list(min_bond_num.keys()))):
                            #print "i",i
                            cutoff=float(min_bond_num[list(min_bond_num.keys())[i]][2])
                            num=min_bond_num[list(min_bond_num.keys())[i]][0]
                            s1=list(min_bond_num.keys())[i]
                            #print s1
                            s2=min_bond_num[list(min_bond_num.keys())[i]][1]
                            s1a=[]
                            s2a=[]
                            for a in range(len(trial_atoms)):
                                if trial_atoms[a].symbol==s1:
                                    s1a.append(a)
                                if trial_atoms[a].symbol==s2:
                                    s2a.append(a)
                            
                            for x in range(len(s1a)):
                                rij=[]
                                for y in range(len(s2a)):
                                    tmp=trial_atoms.get_distance(s1a[x],s2a[y],mic=True)
                                    if tmp <= cutoff and tmp >= min_bond:
                                        rij.append(tmp)
                                bonds.append(rij)
                        for x in range(len(bonds)):
                            #print len(bonds[x])
                            if float(len(bonds[x])) < float(num):
                                test=1
            if test == 1:
                k = k+1
                #print "2  ",k
                if p_t == 'T5' or 'T6':
                    SL = old_sl
                if k >= kmax:
                    print("no more new permutations found")
                    break
                continue
                
        log_file['move_type'].append(p_t)
        log_file['modules'].append(struct)

        ### calculation of energy for new pertubation###
        if len(struct) != SL:
            SL = len(struct)
            print("fault with stacking when using %s, resetting SL to length of current struct obj"%(str(p_t)))
        lowest_energy = energies[(min(energies,key=energies.get))]
        #flag
        #print struct
        new_atoms = perovskite_stack(AMods,BMods,SL,struct,ratt_dist)
        new_atoms,new_energy,converged=run_gulp(atoms=new_atoms,shel=shel,kwds=kwds,opts=opts,lib=lib,gulp_command=gulp_command)
        if converged == False:
            new_energy = 1.0E20
        if graph_out=='T':
            graph.write(str(str(s+p_s).ljust(8)+str('{0:6s}'.format(str("NE = ")))+str('{0:4.8f}'.format(new_energy/len(atoms)).ljust(35))+str(" eV/atom")))
        k = 0
        
        ### insert safe guard here to look at interatomic distances in the newly relaxed structure, reject structure if any distance is equal to or less than 1 angstrom, useful with gulp to prevent structures on the verge of collapsing being accepted as the lowest energy structure. ###
        distances=get_distances(new_atoms=new_atoms)
        if float(min(distances))<= float(min_bond):
            new_energy = 1.0E20
#####################################################################################################
        #if not calc.converged:
        #    new_energy = 1.0E20
        try:
            new_energy=new_energy/len(new_atoms)
        except(TypeError,ValueError):
            new_energy = 1.0E20/len(new_atoms)
        energies[(s+p_s)]=float(new_energy)
        log_file['energy'].append(new_energy)
        diff=float(new_energy-lowest_energy)
        if s >> 1:
            l_conf_old=l_conf
        l_conf=min(energies, key=energies.get)
        if s ==1:
            l_conf_old=l_conf
        outcome = ''
        ### accept new structure if it becomes the lowest energy configuration with chance for MC accept ###
        if l_conf == (s+p_s):
            l_conf_atoms=new_atoms.copy()
            lowest_energy_structure=new_atoms.copy()
            lowest_energy=float(energies[l_conf])
            write(lout,lowest_energy_structure)
            if output == 'full':
                #write("accepted/"+str(s+p_s)+".cif",l_conf_atoms)
                if archive_structures != True:
                    write("accepted/"+str('{0:0=5n}'.format(s+p_s)+".cif"),l_conf_atoms)
            accepted.append(struct)
            outcome = 'A'
            if outcome=='A':
                if new_energy < 1.0e+5:
                    acpt.write(str("\n"+str(s+p_s).ljust(5)+str('{0:4s}'.format(str("NE = ")))+str('{0:4.8f}'.format(new_energy).ljust(35))+str(" eV/atom")))
            r = 0
        if l_conf != (s+p_s):
            p=random.random()
        #### if new energy == old energy, 50% chance to accept new structure as current
            if diff == 0.:
                if p_t == 'T5' or 'T6':
                    chance = 1.
                else:
                    chance = 0.5
                if chance >= p:
                    outcome = 'mcA'
                    l_conf_atoms=new_atoms.copy()
                    lowest_energy=energies[(s+p_s)]
                    #write(lout,l_conf_atoms)
                    accepted.append(struct)
                    if output == 'full':
                        #write("accepted/"+str(s+p_s)+".cif",l_conf_atoms)
                        if archive_structures != True:
                            write("accepted/"+str('{0:0=5n}'.format(s+p_s)+"MCA.cif"),l_conf_atoms)
                        if new_energy < 1.0e+5:
                            acpt.write(str("\n"+str(s+p_s).ljust(5)+str('{0:4s}'.format(str("NE = ")))+str('{0:4.8f}'.format(new_energy).ljust(35))+str(" eV/atom")+str(" MCA")))
                    r = r + 1
                else:    
                    outcome = 'R'
                    r = r+1                
        #### monte-carlo accept test
            elif math.exp(-diff/(float(red_T))) >= p:
                outcome = 'mcA'
                l_conf_atoms=new_atoms.copy()
                lowest_energy=energies[(s+p_s)]
                #write(lout,l_conf_atoms)
                accepted.append(struct)
                if output == 'full':
                    if archive_structures != True:
                        write("accepted/"+str('{0:0=5n}'.format(s+p_s)+"MCA.cif"),l_conf_atoms)
                    #write("accepted/"+str(s+p_s)+"MCA.cif",l_conf_atoms)
                    if new_energy < 1.0e+5:
                        acpt.write(str("\n"+str(s+p_s).ljust(5)+str('{0:4s}'.format(str("NE = ")))+str('{0:4.8f}'.format(new_energy).ljust(35))+str(" eV/atom")+str(" MCA")))

                #if output == 'full':
                #    write("accepted/"+str(s+p_s)+".cif",l_conf_atoms)
                r = r + 1
        
            elif converged == False:
                outcome = 'cR'
                r = r+1
            else:
                outcome = 'R'
                r = r+1


        log_file['mc_outcome'].append(outcome)
        
        if graph_out=='T':
            graph.write(str("\n"+str(s+p_s).ljust(8)+str('{0:6s}'.format(str("NE = ")))+str('{0:4.8f}'.format(new_energy).ljust(35))+str(" eV/atom")))
            
        #print perms
        perms.append(struct)
        SL = len(struct)
        if all_out == 'T':
            #filename = str("steps/step_"+str(s+p_s)+".cif")
            if archive_structures != True:
                write("steps/"+str('{0:0=5n}'.format(s+p_s)+".cif"),new_atoms)
                
            if archive_structures == True:    
                all_structures[s+p_s]=new_atoms
            #write(filename,new_atoms)
        t_s = time.time() - t_i
        if s % 10==0:
            print("")
        #print str(s+p_s), " NE = %6.4e"%(new_energy), " eV/atom, dE = %6.4e"%(diff),"vs.",l_conf_old,"pert =",p_t,"Out:",outcome,"CR:",r,
        print(str(str(s+p_s).ljust(5)+str('{0:4s}'.format(str("NE = ")))+str('{0: 6.4e}'.format(new_energy).ljust(11))+str(" eV/atom dE = {0: 6.4e}".format(diff).ljust(6))+str(" vs.")+str('{0:4n}'.format(l_conf_old).ljust(4))+str(" pert:")+str(p_t)+str(" Out: {0:3s}".format(outcome).ljust(3))+str(" CR: {0:3n}".format(r).ljust(5))), end=' ')

        if dump =='T':
            o.write(str(str(s+p_s).ljust(5)+str('{0:4s}'.format(str("NE = ")))+str('{0: 6.4e}'.format(new_energy).ljust(11))+str(" eV/atom dE = {0: 6.4e}".format(diff).ljust(6))+str(" vs.")+str('{0:4n}'.format(l_conf_old).ljust(4))+str(" pert:")+str(p_t)+str(" Out: {0:3s}".format(outcome).ljust(3))+str(" CR: {0:3n}".format(r).ljust(5))))
        
        ### flag1####
        if len(sl) > 1: #### this doesn't quite make sense at the moment, since sl should stll be a dictionary, so I've messed up somewhere?
#            print "SL =",SL,
            print(str(" SL = {0:3n}".format(len(struct)).ljust(3)), end=' ')
    
            o.write(str(" SL = {0:3n}".format(len(struct)).ljust(3)))
#        print "CPU: %3.1f"%(t_s)
        print(str(" CPU: {0:3.1f}".format(t_s)).ljust(9))
        o.write(str(" CPU: {0:3.1f} \n".format(t_s)).ljust(9))

        o.flush()
        acpt.flush()
        if graph_out == 'T':
            graph.flush()
        sys.stdout.flush()



        ### checking for break criteria/ triggers for delays ###
        if s == smax:
            print("\n*** maximum number of steps reached ***")
        if r >= rmax and (s+p_s) >= smin:
            print("\n*** requested minimum number of consecutive rejections reached ***")
            break
        if p_t in ['T3','T4','T5','T6','T7'] and outcome in ['A','mcA']:
            force = delay 
            f_type='low'
        
        if r >= (rmax *0.9):
            force = 200
            f_type ='high'
        if os.path.isfile("stop.txt"):
            print("soft stop found, aborting run")
            os.remove("stop.txt")
            break
        tmplowest = perms[l_conf].copy()
        eng=[]
        #print energies
        for i in range(len(energies)):
            tmp = [i,energies[i]]
            eng.append(tmp)
        tmpenergies = numpy.array(eng)
        tmpaccepted = numpy.array(accepted)
        numpy.savez("restart",perms=perms,lowest=tmplowest,energies=tmpenergies,r=r,accepted=tmpaccepted,cstack=stack)
        #print(log_file)
        try:
            out_data=pandas.DataFrame(log_file)
        except(ValueError):
            print(log_file)
            sys.exit()
        out_data.to_csv("log_file.csv",index=False)
#        numpy.savez("permitted.npz",permitted=permitted)
        #outfile=open('permitted.dat','wb')
        #pickle.dump(permitted,outfile,protocol=pickle.HIGHEST_PROTOCOL)
        #outfile=open('Sl.dat','wb')
        #pickle.dump(sl,outfile,protocol=pickle.HIGHEST_PROTOCOL)
        s = s+1


    ##### output section ###
    #l_conf=min(energies, key=energies.get)
    write(lout,lowest_energy_structure)
    
    perms = numpy.array(perms)
    print("\nlowest energy configuration is step:",l_conf,"energy = %4.8f eV/atom"%energies[l_conf], ",structures relaxed:",len(perms),"\n")
    o.write(str("\n\nlowest energy configuration is step: {0:n} energy = {1:4.8f} eV/atom, structures relaxed: {2:n}".format(l_conf,energies[l_conf],len(perms))+"\n\n"))
    print("total time :%10.1f"%(time.time()-t_i))



    ###### making restart files: note that energies is converted to an array so that all of the required restart infomation can easily be stored in one .npz file, and is readily converted back to an array at the re-start of a calculation 
    if output == 'full':
        print("writing restart files")
    tmplowest = perms[l_conf].copy()
    eng=[]
    #print energies
    for i in range(len(energies)):
        tmp = [i,energies[i]]
        eng.append(tmp)
    tmpenergies = numpy.array(eng)
    tmpaccepted = numpy.array(accepted)
    if archive_structures == True:
    	 pickle.dump(all_structures,open("all_structures.p",'wb'))
    numpy.savez("restart",perms=perms,lowest=tmplowest,energies=tmpenergies,r=r,accepted=accepted,cstack=stack)
#    numpy.savez("permitted.npz",permitted=permitted)
    #outfile=open('permitted.dat','wb')
    #pickle.dump(permitted,outfile,protocol=pickle.HIGHEST_PROTOCOL)
    #outfile=open('Sl.dat','wb')
    #pickle.dump(sl,outfile,protocol=pickle.HIGHEST_PROTOCOL)
    o.close()
    sys.exit()

