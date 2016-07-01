#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  PS-P2VP_random.py sys
import math
import pylab
import matplotlib.pyplot as plt
import numpy as np
import random

#filename = "config.out"
nA = 8
nB = 8

rab = 10
#nAB = nA+nB
#density = 0.7
number_ions = 0
mpb = 8 #monomer per bead
#sf = 10 #scale factor
bondlength = 1.5
func = 0.5
xl = -5
xh = 6
yl = -5
yh = 4
zl = -5
zh = 4
totm = 0
nchainn = 0
nPolymer = 16

#return the lenth of vector a
def lent(a):
    g = math.sqrt(a.dot(a))
    return g

#return the direction from point b to point a 
def dir(a,b):
    g = (a-b)/lent(a-b)
    return g

#return the cross product of vector a and b
def cross(a,b):
    g = np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[0]-a[1]*b[0]])
    return g

#return the length of the cross product of vector a and b   
def lencross(a,b):
    g = lent(a)*lent(b)*math.sin(math.acos(a.dot(b)/(lent(a)*lent(b))))
    return g
    
#return the unit vector which is vertical to a    
def dirv(a):
    
    if a[2] !=0 :
        g = np.array([0,1,0])
        g[2] = -a[1]/a[2]
    else:
        g = np.array([0,1,0])
        g[0] = -a[1]/a[0] 
    g = g/lent(g)
    return g
    
#return the unit vector which is vertical to both a and b
def dirvd(a,b):
    c = np.asarray([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])
    c = c/lent(c)
    return c

#return the distance between point a and b
def lenp(a,b):
    g = np.array([b[0]-a[0],b[1]-a[1],b[2]-a[2]])
    m = lent(g)
    return m
    
#construct PS
def getposh_ali1(cali1, vec):
    x = cali1[0] - bondlength*vec[0]
    y = cali1[1] - bondlength*vec[1]
    z = cali1[2] - bondlength*vec[2]
    pos = np.asarray([x,y,z])
    return pos  
    
def getposh_ali2(cali1, vecv):
    x = cali1[0] + bondlength*vecv[0]
    y = cali1[1] + bondlength*vecv[1]
    z = cali1[2] + bondlength*vecv[2]
    pos = np.asarray([x,y,z])
    return pos  

def getposh_ali3(cali1, vecv):
    x = cali1[0] - bondlength*vecv[0]
    y = cali1[1] - bondlength*vecv[1]
    z = cali1[2] - bondlength*vecv[2]
    pos = np.asarray([x,y,z])
    return pos  

def getposc_ali2(cali1, vecm):
    x = cali1[0] + bondlength*vecm[0]
    y = cali1[1] + bondlength*vecm[1]
    z = cali1[2] + bondlength*vecm[2]
    pos = np.asarray([x,y,z])
    return pos  
    
def getposh_ali4(cali2,vecv):
    if random.random()>0.5:
        x = cali2[0] + bondlength*vecv[0]
        y = cali2[1] + bondlength*vecv[1]
        z = cali2[2] + bondlength*vecv[2]
    else:
        x = cali2[0] - bondlength*vecv[0]
        y = cali2[1] - bondlength*vecv[1]
        z = cali2[2] - bondlength*vecv[2]
    pos = np.asarray([x,y,z])
    return pos

def getposc_aro6(cali2,vecd):
    x = cali2[0] + bondlength*vecd[0]
    y = cali2[1] + bondlength*vecd[1]
    z = cali2[2] + bondlength*vecd[2]
    pos = np.asarray([x,y,z])
    return pos

def getposc_aro5(caro6, vecmr):
    x = caro6[0] + bondlength*vecmr[0]
    y = caro6[1] + bondlength*vecmr[1]
    z = caro6[2] + bondlength*vecmr[2]
    pos = np.asarray([x,y,z])
    return pos
    
def getposh_aro5(caro5,vecm):
    x = caro5[0] - bondlength*vecm[0]
    y = caro5[1] - bondlength*vecm[1]
    z = caro5[2] - bondlength*vecm[2]
    pos = np.asarray([x,y,z])
    return pos

def getposc_aro4(caro5, vecd):
    x = caro5[0] + bondlength*vecd[0]
    y = caro5[1] + bondlength*vecd[1]
    z = caro5[2] + bondlength*vecd[2]
    pos = np.asarray([x,y,z])
    return pos
    
def getposh_aro4(caro4, vecmr):
    x = caro4[0] + bondlength*vecmr[0]
    y = caro4[1] + bondlength*vecmr[1]
    z = caro4[2] + bondlength*vecmr[2]
    pos = np.asarray([x,y,z])
    return pos

def getposc_aro3(caro6,vecd):
    x = caro6[0] + 2*bondlength*vecd[0]
    y = caro6[1] + 2*bondlength*vecd[1]
    z = caro6[2] + 2*bondlength*vecd[2]
    pos = np.asarray([x,y,z])
    return pos

def getposh_aro3(caro3, vecd):
    x = caro3[0] + bondlength*vecd[0]
    y = caro3[1] + bondlength*vecd[1]
    z = caro3[2] + bondlength*vecd[2]
    pos = np.asarray([x,y,z])
    return pos

def getposc_aro1(caro6,vecm):
    x = caro6[0] + bondlength*vecm[0]
    y = caro6[1] + bondlength*vecm[1]
    z = caro6[2] + bondlength*vecm[2]   
    pos = np.asarray([x,y,z])
    return pos

def getposh_aro1(caro1, vecmr):
    x = caro1[0] - bondlength*vecmr[0]
    y = caro1[1] - bondlength*vecmr[1]
    z = caro1[2] - bondlength*vecmr[2]
    pos = np.asarray([x,y,z])
    return pos
    
def getposc_aro2(caro1,vecd):
    x = caro1[0] + bondlength*vecd[0]
    y = caro1[1] + bondlength*vecd[1]
    z = caro1[2] + bondlength*vecd[2]
    pos = np.asarray([x,y,z])
    return pos

def getposh_aro2(caro2,vecm):
    x = caro2[0] + bondlength*vecm[0]
    y = caro2[1] + bondlength*vecm[1]
    z = caro2[2] + bondlength*vecm[2]
    pos = np.asarray([x,y,z])
    return pos
    
def write_data_atom(index_atom, atom):
    #x_coord = coord[0]
    #y_coord = coord[1]
    #z_coord = coord[2]
    #print (atom)
    index_chain = atom[0]
    type_atom = atom[1]
    charge = atom[2]
    x_coord = atom[3][0]
    y_coord = atom[3][1]
    z_coord = atom[3][2]
    
    index_chain +=1
    #print#print len(name)
    with open("PS_P2VP_random.data", "a") as out_file:
        #out_file.write(str(index_atom) + ' ' + str(atom).strip('[]').strip(',') + '\n')
        out_file.write(str(index_atom)+ ' ' + str(index_chain)+ ' ' + str(type_atom) + ' ' + str(charge) + ' ' + 
        str("{0:.3f}".format(x_coord)).strip('[]') + ' ' + str("{0:.3f}".format(y_coord)).strip('[]') + ' ' +
        str("{0:.3f}".format(z_coord)).strip('[]') + '\n')#+ ' ' + name + ' ' + monomer)

def write_data_bond(index_bond, type_bond):#, atom1, atom2):
    with open("PS_P2VP_random.data", "a") as out_file:
        out_file.write(str(index_bond) + ' ' + str(type_bond).strip('[]').strip(',')+'\n')# + ' ' + str(atom1) + ' ' + str(atom2) + '\n')

def write_data_angle(index_angle, type_angle):#, atom1, atom2, atom3):
    with open("PS_P2VP_random.data", "a") as out_file:
        out_file.write(str(index_angle) + ' ' + str(type_angle).strip('[]').strip(',') + '\n')# + ' ' + str(atom1) + ' ' + str(atom2) + ' ' + str(atom3) + '\n')

def write_data_dihedral(index_dihedral, type_dihedral):#, atom1, atom2, atom3, atom4):
    with open("PS_P2VP_random.data", "a") as out_file:
        out_file.write(str(index_dihedral) + ' ' + str(type_dihedral).strip('[]').strip(',') + '\n')# + ' ' + str(atom1) + ' ' + str(atom2) + ' ' + str(atom3) + ' ' + str(atom4) + '\n')

def write_data_improper(index_improper, type_improper):#, atom1, atom2, atom3, atom4):
    with open("PS_P2VP_random.data", "a") as out_file:
        out_file.write(str(index_improper) + ' ' + str(type_improper).strip('[]').strip(',') + '\n')# + ' ' + str(atom1) + ' ' + str(atom2) + ' ' + str(atom3) + ' ' + str(atom4) + '\n')   
                            

d = open('sconfig.xyz', 'w')
#m = int(math.sqrt(nPolymer))

for i in range(nPolymer):
	for j in range(nA):
		if j == 0:
			x = random.random()*rab
			y = random.random()*rab
			z = random.random()*rab
		else:
			theta = random.random()*2*math.pi
			if random.random()>0.5:
				flag = 1
			else:
				flag = -1
			dz = random.random()*flag
			
			dx = math.sqrt(1-dz*dz)*math.cos(theta)
			dy = math.sqrt(1-dz*dz)*math.sin(theta)
			dz = dz
			x += dx
			y += dy
			z += dz
		d.write('{0:8d}{1:15.5f}{2:15.5f}{3:15.5f}\n'.format(int(0), x, y, z))	
	for j in range(nB):
		theta = random.random()*2*math.pi
		if random.random()>0.5:
			flag = 1
		else:
			flag = -1
		dz = random.random()*flag
		
		dx = math.sqrt(1-dz*dz)*math.cos(theta)
		dy = math.sqrt(1-dz*dz)*math.sin(theta)
		dz = dz
		x += dx
		y += dy
		z += dz
		d.write('{0:8d}{1:15.5f}{2:15.5f}{3:15.5f}\n'.format(int(1), x, y, z))

          
d.close()
t0 = []
t1 = []

d = open('sconfig.xyz','r')
lines=d.readlines()
for i in range(len(lines)):
    lines[i] = lines[i].split()
    lines[i] = np.array(map(float, lines[i]))
    if(lines[i][0] == 0):
        t0.append([lines[i][1],lines[i][2],lines[i][3]])
    else:
		t1.append([lines[i][1],lines[i][2],lines[i][3]])


k = 0 
bond = 0      
for i in range(1,len(t0)):
    if(i%nA != 0):
        bond += lenp(t0[i-1],t0[i])
        k += 1
bond /= k
print "bondlength for coarse grain bead  ", bond

#rescale the bead position to make the bondlength = bondlength
#unit = angstroms
sf = bondlength * 1.732 * mpb/bond
for i in range(len(t0)):
    t0[i] = [t0[i][0]*sf, t0[i][1]*sf, t0[i][2]*sf]
for i in range(len(t1)):
    t1[i] = [t1[i][0]*sf, t1[i][1]*sf, t1[i][2]*sf]
    
k = 0 
bond = 0     
for i in range(1,len(t0)):
    if(i%nA != 0):
        bond += lenp(t0[i-1],t0[i])
        k += 1
bond /= k
print "bondlength for scaled coarse grain bead  ", bond
#bondlength between monomers
bondm = 1.732*bondlength


#back mapping N polymer chains with nA A beads and nB B beads, one bead 
#will be back mapped to nM monomers

num_chains = nPolymer
#num_int_monomers_P2VP = nA -1
#num_int_monomers_PMMA = nB - 1
index_monomer = 1


index_atom = 1
index_chain = 1

atom = []
bonding = []
angle = []
dihedral = []
improper = []
y_count = 0
z_count = 0
number_ions = 0 
k = 0
j = 0
s = 0   

for i in range(num_chains):
    #print "chain",i
    backbone_dihedral = []
    for k in range(nA):
        if(k == nA-1):
            mpb = int(lenp(t1[i*nB],t0[i*nA+k])/bondm)
        else:
            mpb = int(lenp(t0[i*nA+k+1],t0[i*nA+k])/bondm)
            #print lenp(t0[i*nA+k+1],t0[i*nA+k])
            #print mpb
        totm += mpb
        for j in range(mpb):
            if(k == 0 and j == 0):
                #print " starting monomer", index_atom
                #print len(atom),"len(atom)"
               # print " k==0,j==0", k,j
                #Initial monomer of PS
                vec = np.asarray([t0[i*nA+k+1][0]-t0[i*nA+k][0],t0[i*nA+k+1][1]-t0[i*nA+k][1],t0[i*nA+k+1][2]-t0[i*nA+k][2]])
                #print k, j
               # print "vec", vec
                lentv = lent(vec)
                #the vector which is the direction of two bead
                vec = vec/lent(vec)
                #a unit vector which is vertical t vec
               # print vec
                vecv = dirv(vec)
        
                vecd = dirvd(vec, vecv)
                #print vecd
                vecm = 1.732*vec + vecd
                vecm = vecm/lent(vecm)
                vecmr = vecd - 1.732*vec
                vecmr = vecmr/lent(vecmr)
                vecmv = 1.732*vecv +vecd
                vecmv = vecmv/lent(vecmv)
                vecmrv = vecd - 1.732*vecv
                vecmrv = vecmrv/lent(vecmrv)
                dx = vec[0]*lentv/mpb
                dy = vec[1]*lentv/mpb
                dz = vec[2]*lentv/mpb           
                x = t0[i*nB+k][0]
                y = t0[i*nB+k][1]
                z = t0[i*nB+k][2]
                
                posc_ali1 = np.asarray([x,y,z])
                C_ali1 = np.asarray([x,y,z, index_atom])
                atom.append([i, 1, 0.0, C_ali1])
                index_atom +=1
                            
                posh_ali1 = getposh_ali1(posc_ali1,vec)
                H_ali1 = np.asarray([posh_ali1[0], posh_ali1[1], posh_ali1[2], index_atom])
                atom.append([i, 2, 0.0, H_ali1])
                index_atom +=1
                
                posh_ali2 = getposh_ali2(posc_ali1, vecv)
                H_ali2 = np.asarray([posh_ali2[0], posh_ali2[1], posh_ali2[2], index_atom])
                atom.append([i, 2, 0.0, H_ali2])
                index_atom +=1
                
                posh_ali3 = getposh_ali3(posc_ali1, vecv)
                H_ali3 = np.asarray([posh_ali3[0], posh_ali3[1], posh_ali3[2], index_atom])
                atom.append([i, 2, 0.0, H_ali3])
                index_atom +=1
                
                posc_ali2 = getposc_ali2(posc_ali1, vecm)
                C_ali2 = np.asarray([posc_ali2[0], posc_ali2[1], posc_ali2[2],index_atom])
                atom.append([i, 1, 0.0, C_ali2])
                index_atom +=1
                
                posh_ali4 = getposh_ali4(posc_ali2,vecv)
                H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                atom.append([i, 2, 0.0, H_ali4])
                index_atom +=1
                
                posc_aro6 = getposc_aro6(posc_ali2,vecd)
                C_aro6 = np.asarray([posc_aro6[0], posc_aro6[1], posc_aro6[2], index_atom])
                atom.append([i, 3, 0.0, C_aro6])
                index_atom +=1
                
                posc_aro5 = getposc_aro5(posc_aro6,vecmrv)
                C_aro5 = np.asarray([posc_aro5[0], posc_aro5[1], posc_aro5[2], index_atom])
                atom.append([i, 4, -0.115, C_aro5])
                index_atom +=1
                
                posh_aro5 = getposh_aro5(posc_aro5, vecmv)
                H_aro5 = np.asarray([posh_aro5[0], posh_aro5[1], posh_aro5[2], index_atom])
                atom.append([i, 5, 0.115, H_aro5])
                index_atom +=1
                
                posc_aro4 = getposc_aro4(posc_aro5,vecd)
                C_aro4 = np.asarray([posc_aro4[0], posc_aro4[1], posc_aro4[2], index_atom])
                atom.append([i, 4, -0.115, C_aro4])
                index_atom +=1
                
                posh_aro4 = getposh_aro4(posc_aro4,vecmrv)
                H_aro4 = np.asarray([posh_aro4[0], posh_aro4[1], posh_aro4[2], index_atom])
                atom.append([i, 5, 0.115, H_aro4])
                index_atom +=1
                
                posc_aro3 = getposc_aro3(posc_aro6,vecd)
                C_aro3 = np.asarray([posc_aro3[0], posc_aro3[1], posc_aro3[2], index_atom])
                atom.append([i, 4, -0.115, C_aro3])
                index_atom +=1
                
                posh_aro3 = getposh_aro3(posc_aro3,vecd)
                H_aro3 = np.asarray([posh_aro3[0], posh_aro3[1], posh_aro3[2], index_atom])
                atom.append([i, 5, 0.115, H_aro3])
                index_atom +=1
                
                posc_aro1 = getposc_aro1(posc_aro6, vecmv)
                C_aro1 = np.asarray([posc_aro1[0], posc_aro1[1], posc_aro1[2], index_atom])
                atom.append([i, 4, -0.115, C_aro1])
                index_atom +=1
                
                posh_aro1 = getposh_aro1(posc_aro1, vecmrv)
                H_aro1 = np.asarray([posh_aro1[0], posh_aro1[1], posh_aro1[2], index_atom])
                atom.append([i, 5, 0.115, H_aro1])
                index_atom +=1
                
                posc_aro2 = getposc_aro2(posc_aro1, vecd)
                C_aro2 = np.asarray([posc_aro2[0], posc_aro2[1], posc_aro2[2], index_atom])
                atom.append([i, 4, -0.115, C_aro2])
                index_atom +=1
                
                posh_aro2 = getposh_aro2(posc_aro2, vecm)
                H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                atom.append([i, 5, 0.115, H_aro2])
                index_atom +=1
                
                #Storing bonding information
                bonding.append([4, int(C_ali1[3]), int(H_ali1[3])])
                bonding.append([4, int(C_ali1[3]), int(H_ali2[3])])
                bonding.append([4, int(C_ali1[3]), int(H_ali3[3])])            
                bonding.append([1, int(C_ali1[3]), int(C_ali2[3])])
                bonding.append([4, int(C_ali2[3]), int(H_ali4[3])])
                bonding.append([2, int(C_ali2[3]), int(C_aro6[3])])
                bonding.append([3, int(C_aro5[3]), int(C_aro6[3])])
                bonding.append([5, int(C_aro5[3]), int(H_aro5[3])])
                bonding.append([3, int(C_aro4[3]), int(C_aro5[3])])
                bonding.append([5, int(C_aro4[3]), int(H_aro4[3])])
                bonding.append([3, int(C_aro3[3]), int(C_aro4[3])])
                bonding.append([5, int(C_aro3[3]), int(H_aro3[3])])
                bonding.append([3, int(C_aro1[3]), int(C_aro6[3])])
                bonding.append([5, int(C_aro1[3]), int(H_aro1[3])])
                bonding.append([3, int(C_aro2[3]), int(C_aro1[3])])
                bonding.append([3, int(C_aro2[3]), int(C_aro3[3])])
                bonding.append([5, int(C_aro2[3]), int(H_aro2[3])])
                bonding.append([1, int(C_ali2[3]), int(index_atom)]) ##bonding to next monomer
                
                ##Storing angle information of monomer
                angle.append([2, int(H_ali1[3]), int(C_ali1[3]), int(H_ali2[3])])
                angle.append([2, int(H_ali1[3]), int(C_ali1[3]), int(H_ali3[3])])
                angle.append([2, int(H_ali3[3]), int(C_ali1[3]), int(H_ali2[3])])
                angle.append([1, int(H_ali1[3]), int(C_ali1[3]), int(C_ali2[3])])
                angle.append([1, int(H_ali2[3]), int(C_ali1[3]), int(C_ali2[3])])
                angle.append([1, int(H_ali3[3]), int(C_ali1[3]), int(C_ali2[3])])
                angle.append([1, int(H_ali4[3]), int(C_ali2[3]), int(C_ali1[3])])
                angle.append([1, int(H_ali4[3]), int(C_ali2[3]), int(index_atom)])
                angle.append([1, int(index_atom+1), int(index_atom), int(C_ali2[3])])
                angle.append([1, int(index_atom+2), int(index_atom), int(C_ali2[3])])
                angle.append([1, int(C_aro6[3]), int(C_ali2[3]), int(H_ali4[3])])
                angle.append([3, int(C_ali1[3]), int(C_ali2[3]), int(index_atom)])
                angle.append([3, int(C_ali2[3]), int(index_atom), int(index_atom+3)])
                angle.append([3, int(C_ali1[3]), int(C_ali2[3]), int(C_aro6[3])])
                angle.append([3, int(index_atom), int(C_ali2[3]), int(C_aro6[3])])
                angle.append([4, int(C_ali2[3]), int(C_aro6[3]), int(C_aro5[3])])
                angle.append([4, int(C_ali2[3]), int(C_aro6[3]), int(C_aro1[3])])
                angle.append([4, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                angle.append([4, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                angle.append([4, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                angle.append([4, int(C_aro3[3]), int(C_aro2[3]), int(C_aro1[3])])
                angle.append([4, int(C_aro2[3]), int(C_aro1[3]), int(C_aro6[3])])
                angle.append([4, int(C_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                angle.append([5, int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                angle.append([5, int(C_aro4[3]), int(C_aro5[3]), int(H_aro5[3])])
                angle.append([5, int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                angle.append([5, int(C_aro3[3]), int(C_aro4[3]), int(H_aro4[3])])
                angle.append([5, int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                angle.append([5, int(C_aro2[3]), int(C_aro3[3]), int(H_aro3[3])])
                angle.append([5, int(C_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                angle.append([5, int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                angle.append([5, int(C_aro6[3]), int(C_aro1[3]), int(H_aro1[3])])
                angle.append([5, int(C_aro2[3]), int(C_aro1[3]), int(H_aro1[3])])
            
        ##Storing dihedral information of monomer
                #dihedral.append([1, int(C_ali1[3]), int(C_ali2[3]), int(index_atom), int(index_atom+3)])
                backbone_dihedral.append(int(C_ali1[3]))
                backbone_dihedral.append(int(C_ali2[3]))
                #dihedral.append([1, int(C_ali2[3]), int(index_atom), int(index_atom+3), int(index_atom+16)])
                dihedral.append([1, int(H_ali1[3]), int(C_ali1[3]), int(C_ali2[3]), int(index_atom)])
                dihedral.append([1, int(H_ali2[3]), int(C_ali1[3]), int(C_ali2[3]), int(index_atom)])
                dihedral.append([1, int(H_ali3[3]), int(C_ali1[3]), int(C_ali2[3]), int(index_atom)])
            
        ##Storing improper information of monomer
                improper.append([1, int(C_aro6[3]), int(C_aro5[3]), int(C_aro1[3]), int(C_ali2[3])])
                improper.append([1, int(C_aro5[3]), int(C_aro4[3]), int(C_aro6[3]), int(H_aro5[3])])
                improper.append([1, int(C_aro4[3]), int(C_aro3[3]), int(C_aro5[3]), int(H_aro4[3])])
                improper.append([1, int(C_aro3[3]), int(C_aro2[3]), int(C_aro4[3]), int(H_aro3[3])])
                improper.append([1, int(C_aro2[3]), int(C_aro3[3]), int(C_aro1[3]), int(H_aro2[3])])
                improper.append([1, int(C_aro1[3]), int(C_aro2[3]), int(C_aro6[3]), int(H_aro1[3])])
                improper.append([1, int(C_aro1[3]), int(C_aro2[3]), int(C_aro3[3]), int(C_aro4[3])])
                improper.append([1, int(C_aro2[3]), int(C_aro3[3]), int(C_aro4[3]), int(C_aro5[3])])
                improper.append([1, int(C_aro3[3]), int(C_aro4[3]), int(C_aro5[3]), int(C_aro6[3])])
                improper.append([1, int(C_aro4[3]), int(C_aro5[3]), int(C_aro6[3]), int(C_aro1[3])])
                improper.append([1, int(C_aro5[3]), int(C_aro6[3]), int(C_aro1[3]), int(C_aro2[3])])
                improper.append([1, int(C_aro6[3]), int(C_aro1[3]), int(C_aro2[3]), int(C_aro3[3])])
                
            else:
              #  print"k and j at least one is not zero", k, j
                if (j==0):
                    if(k == nA-1):
                        vec = np.asarray([t1[nB*i][0]-t0[i*nA+k][0],t1[nB*i][1]-t0[i*nA+k][1],t1[nB*i][2]-t0[i*nA+k][2]])
                        #s = s+nB
                    else:   
                        vec = np.asarray([t0[i*nA+k+1][0]-t0[i*nA+k][0],t0[i*nA+k+1][1]-t0[i*nA+k][1],t0[i*nA+k+1][2]-t0[i*nA+k][2]])
                    lentv = lent(vec)
                    #the vector which is the direction of two bead
                    vec = vec/lent(vec)
                    #a unit vector which is vertical t vec
                   # print vec
                    vecv = dirv(vec)
                    vecd = dirvd(vec, vecv)
                    vecm = 1.732*vec + vecd
                    vecm = vecm/lent(vecm)
                    vecmr = vecd - 1.732*vec
                    vecmr = vecmr/lent(vecmr)
                    vecmv = 1.732*vecv +vecd
                    vecmv = vecmv/lent(vecmv)
                    vecmrv = vecd - 1.732*vecv
                    vecmrv = vecmrv/lent(vecmrv)
                    dx = vec[0]*lentv/mpb
                    dy = vec[1]*lentv/mpb
                    dz = vec[2]*lentv/mpb 
                   # print "bondlenth", math.sqrt(dx*dx+dy*dy+dz*dz), "should be 2.5"
                    x = t0[i*nA+k][0]
                    y = t0[i*nA+k][1]
                    z = t0[i*nA+k][2] 
                else:         
                    x = x + dx
                    y = y + dy
                    z = z + dz
                
                
                posc_ali1 = np.asarray([x,y,z])
                C_ali3 = np.asarray([x,y,z, index_atom])
                atom.append([i, 1, 0.0, C_ali3])
                index_atom +=1
                
                C_ali1 = C_ali3
                            
                posh_ali4 = getposh_ali2(posc_ali1, vecv)
                H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                atom.append([i, 2, 0.0, H_ali4])
                index_atom +=1
                
                posh_ali5 = getposh_ali3(posc_ali1, vecv)
                H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                atom.append([i, 2, 0.0, H_ali5])
                index_atom +=1
                
                posc_ali2 = getposc_ali2(posc_ali1, vecm)
                C_ali4 = np.asarray([posc_ali2[0], posc_ali2[1], posc_ali2[2],index_atom])
                atom.append([i, 1, 0.0, C_ali4])
                index_atom +=1
                
                C_ali2 = C_ali4
                
                posh_ali4 = getposh_ali4(posc_ali2,vecv)
                H_ali3 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                atom.append([i, 2, 0.0, H_ali3])
                index_atom +=1
                
                posc_aro6 = getposc_aro6(posc_ali2,vecd)
                C_aro6 = np.asarray([posc_aro6[0], posc_aro6[1], posc_aro6[2], index_atom])
                atom.append([i, 3, 0.0, C_aro6])
                index_atom +=1
                
                posc_aro5 = getposc_aro5(posc_aro6,vecmrv)
                C_aro5 = np.asarray([posc_aro5[0], posc_aro5[1], posc_aro5[2], index_atom])
                atom.append([i, 4, -0.115, C_aro5])
                index_atom +=1
                
                posh_aro5 = getposh_aro5(posc_aro5, vecmv)
                H_aro5 = np.asarray([posh_aro5[0], posh_aro5[1], posh_aro5[2], index_atom])
                atom.append([i, 5, 0.115, H_aro5])
                index_atom +=1
                
                posc_aro4 = getposc_aro4(posc_aro5,vecd)
                C_aro4 = np.asarray([posc_aro4[0], posc_aro4[1], posc_aro4[2], index_atom])
                atom.append([i, 4, -0.115, C_aro4])
                index_atom +=1
                
                posh_aro4 = getposh_aro4(posc_aro4,vecmrv)
                H_aro4 = np.asarray([posh_aro4[0], posh_aro4[1], posh_aro4[2], index_atom])
                atom.append([i, 5, 0.115, H_aro4])
                index_atom +=1
                
                posc_aro3 = getposc_aro3(posc_aro6,vecd)
                C_aro3 = np.asarray([posc_aro3[0], posc_aro3[1], posc_aro3[2], index_atom])
                atom.append([i, 4, -0.115, C_aro3])
                index_atom +=1
                
                posh_aro3 = getposh_aro3(posc_aro3,vecd)
                H_aro3 = np.asarray([posh_aro3[0], posh_aro3[1], posh_aro3[2], index_atom])
                atom.append([i, 5, 0.115, H_aro3])
                index_atom +=1
                
                posc_aro1 = getposc_aro1(posc_aro6, vecmv)
                C_aro1 = np.asarray([posc_aro1[0], posc_aro1[1], posc_aro1[2], index_atom])
                atom.append([i, 4, -0.115, C_aro1])
                index_atom +=1
                
                posh_aro1 = getposh_aro1(posc_aro1, vecmrv)
                H_aro1 = np.asarray([posh_aro1[0], posh_aro1[1], posh_aro1[2], index_atom])
                atom.append([i, 5, 0.115, H_aro1])
                index_atom +=1
                
                posc_aro2 = getposc_aro2(posc_aro1, vecd)
                C_aro2 = np.asarray([posc_aro2[0], posc_aro2[1], posc_aro2[2], index_atom])
                atom.append([i, 4, -0.115, C_aro2])
                index_atom +=1
                
                posh_aro2 = getposh_aro2(posc_aro2, vecmv)
                H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                atom.append([i, 5, 0.115, H_aro2])
                index_atom +=1
                
                
                ### Bonding information for monomer
                bonding.append([4, int(C_ali3[3]), int(H_ali4[3])])
                bonding.append([4, int(C_ali3[3]), int(H_ali5[3])])
                bonding.append([1, int(C_ali3[3]), int(C_ali4[3])])
                bonding.append([4, int(C_ali4[3]), int(H_ali3[3])])
                bonding.append([2, int(C_ali4[3]), int(C_aro6[3])])
                bonding.append([3, int(C_aro5[3]), int(C_aro6[3])])
                bonding.append([5, int(C_aro5[3]), int(H_aro5[3])])
                bonding.append([3, int(C_aro4[3]), int(C_aro5[3])])
                bonding.append([5, int(C_aro4[3]), int(H_aro4[3])])
                bonding.append([3, int(C_aro3[3]), int(C_aro4[3])])
                bonding.append([5, int(C_aro3[3]), int(H_aro3[3])])
                bonding.append([3, int(C_aro1[3]), int(C_aro6[3])])
                bonding.append([5, int(C_aro1[3]), int(H_aro1[3])])
                bonding.append([3, int(C_aro2[3]), int(C_aro1[3])])
                bonding.append([3, int(C_aro2[3]), int(C_aro3[3])])
                bonding.append([5, int(C_aro2[3]), int(H_aro2[3])])
                bonding.append([1, int(C_ali4[3]), int(index_atom)])#bonding to next monomer
                
                
        ##Storing angle information of monomer
                angle.append([2, int(H_ali4[3]), int(C_ali3[3]), int(H_ali5[3])])
                angle.append([1, int(H_ali5[3]), int(C_ali3[3]), int(C_ali4[3])])
                angle.append([1, int(H_ali3[3]), int(C_ali4[3]), int(C_ali3[3])])
                angle.append([1, int(H_ali4[3]), int(C_ali3[3]), int(C_ali4[3])])
                angle.append([1, int(H_ali3[3]), int(C_ali4[3]), int(index_atom)])
                angle.append([1, int(index_atom+1), int(index_atom), int(C_ali4[3])])
                angle.append([1, int(index_atom+2), int(index_atom), int(C_ali4[3])])
                angle.append([1, int(C_aro6[3]), int(C_ali4[3]), int(H_ali3[3])])
                angle.append([3, int(C_ali3[3]), int(C_ali4[3]), int(index_atom)])
                angle.append([3, int(C_ali4[3]), int(index_atom), int(index_atom+3)])
                angle.append([3, int(C_ali3[3]), int(C_ali4[3]), int(C_aro6[3])])
                angle.append([3, int(index_atom), int(C_ali4[3]), int(C_aro6[3])])
                angle.append([4, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3])])
                angle.append([4, int(C_ali4[3]), int(C_aro6[3]), int(C_aro1[3])])
                angle.append([4, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                angle.append([4, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                angle.append([4, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                angle.append([4, int(C_aro3[3]), int(C_aro2[3]), int(C_aro1[3])])
                angle.append([4, int(C_aro2[3]), int(C_aro1[3]), int(C_aro6[3])])
                angle.append([4, int(C_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                angle.append([5, int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                angle.append([5, int(C_aro4[3]), int(C_aro5[3]), int(H_aro5[3])])
                angle.append([5, int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                angle.append([5, int(C_aro3[3]), int(C_aro4[3]), int(H_aro4[3])])
                angle.append([5, int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                angle.append([5, int(C_aro2[3]), int(C_aro3[3]), int(H_aro3[3])])
                angle.append([5, int(C_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                angle.append([5, int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                angle.append([5, int(C_aro6[3]), int(C_aro1[3]), int(H_aro1[3])])
                angle.append([5, int(C_aro2[3]), int(C_aro1[3]), int(H_aro1[3])])
                

                backbone_dihedral.append(int(C_ali3[3]))
                backbone_dihedral.append(int(C_ali4[3]))
            
            ##Storing improper information of monomer
                improper.append([1, int(C_aro6[3]), int(C_aro5[3]), int(C_aro1[3]), int(C_ali2[3])])
                improper.append([1, int(C_aro5[3]), int(C_aro4[3]), int(C_aro6[3]), int(H_aro5[3])])
                improper.append([1, int(C_aro4[3]), int(C_aro3[3]), int(C_aro5[3]), int(H_aro4[3])])
                improper.append([1, int(C_aro3[3]), int(C_aro2[3]), int(C_aro4[3]), int(H_aro3[3])])
                improper.append([1, int(C_aro2[3]), int(C_aro3[3]), int(C_aro1[3]), int(H_aro2[3])])
                improper.append([1, int(C_aro1[3]), int(C_aro2[3]), int(C_aro6[3]), int(H_aro1[3])])
                improper.append([1, int(C_aro1[3]), int(C_aro2[3]), int(C_aro3[3]), int(C_aro4[3])])
                improper.append([1, int(C_aro2[3]), int(C_aro3[3]), int(C_aro4[3]), int(C_aro5[3])])
                improper.append([1, int(C_aro3[3]), int(C_aro4[3]), int(C_aro5[3]), int(C_aro6[3])])
                improper.append([1, int(C_aro4[3]), int(C_aro5[3]), int(C_aro6[3]), int(C_aro1[3])])
                improper.append([1, int(C_aro5[3]), int(C_aro6[3]), int(C_aro1[3]), int(C_aro2[3])])
                improper.append([1, int(C_aro6[3]), int(C_aro1[3]), int(C_aro2[3]), int(C_aro3[3])])
                
    for k in range(nB):
        #print "k", k
        #print "i",i
        if(k == nB-1):
            mpb = int(lenp(t1[i*nB+k-1],t1[i*nB+k])/bondm)
        else:
            mpb = int(lenp(t1[i*nB+k+1],t1[i*nB+k])/bondm)
        totm += mpb
        for j in range(mpb):            
            if(j == mpb-1 and k == nB-1): 
               # print nchainn
                #nchainn +=1
                #print index_atom,"after chain"
                x = x + dx
                y = y + dy
                z = z + dz
                funct = random.random()
                if (funct > func):
                    #print "no func"
                    
                    posc_ali1 = np.asarray([x,y,z])
                    C_ali3 = np.asarray([x,y,z, index_atom])
                    atom.append([i, 6, 0.0, C_ali3])
                    index_atom +=1
                    
                    C_ali1 = C_ali3
                                
                    posh_ali4 = getposh_ali2(posc_ali1, vecv)
                    H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali4])
                    index_atom +=1
                    
                    posh_ali5 = getposh_ali3(posc_ali1, vecv)
                    H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali5])
                    index_atom +=1
                    
                    #print index_atom,"last monomer"
                    with open("lastmonomer.txt",'a') as outa:
                        outa.write(str(index_atom)+'\n')
                    
                    posc_ali2 = getposc_ali2(posc_ali1, vecm)
                    C_ali4 = np.asarray([posc_ali2[0], posc_ali2[1], posc_ali2[2],index_atom])
                    atom.append([i, 6, 0.0, C_ali4])
                    index_atom +=1
                    
                    C_ali2 = C_ali4
                    
                    posh_ali4 = getposh_ali4(posc_ali2,vecv)
                    H_ali3 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali3])
                    index_atom +=1
                    
                    posh_ali2 = getposh_ali1(posc_ali2,-vec)
                    H_ali2 = np.asarray([posh_ali2[0], posh_ali2[1], posh_ali2[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali2])
                    index_atom +=1
                    
                    posc_aro6 = getposc_aro6(posc_ali2,vecd)
                    C_aro6 = np.asarray([posc_aro6[0], posc_aro6[1], posc_aro6[2], index_atom])
                    atom.append([i, 8, 0.4742, C_aro6])
                    index_atom +=1
                    
                    posc_aro5 = getposc_aro5(posc_aro6,vecmrv)
                    C_aro5 = np.asarray([posc_aro5[0], posc_aro5[1], posc_aro5[2], index_atom])
                    atom.append([i, 8, -0.4458, C_aro5])
                    index_atom +=1
                    
                    posh_aro5 = getposh_aro5(posc_aro5, vecmv)
                    H_aro5 = np.asarray([posh_aro5[0], posh_aro5[1], posh_aro5[2], index_atom])
                    atom.append([i, 9, 0.1562, H_aro5])
                    index_atom +=1
                    
                    posc_aro4 = getposc_aro4(posc_aro5,vecd)
                    C_aro4 = np.asarray([posc_aro4[0], posc_aro4[1], posc_aro4[2], index_atom])
                    atom.append([i, 8, 0.2282, C_aro4])
                    index_atom +=1
                    
                    posh_aro4 = getposh_aro4(posc_aro4,vecmrv)
                    H_aro4 = np.asarray([posh_aro4[0], posh_aro4[1], posh_aro4[2], index_atom])
                    atom.append([i, 9, 0.0662, H_aro4])
                    index_atom +=1
                    
                    posc_aro3 = getposc_aro3(posc_aro6,vecd)
                    C_aro3 = np.asarray([posc_aro3[0], posc_aro3[1], posc_aro3[2], index_atom])
                    atom.append([i, 8, -0.4458, C_aro3])
                    index_atom +=1
                    
                    posh_aro3 = getposh_aro3(posc_aro3,vecd)
                    H_aro3 = np.asarray([posh_aro3[0], posh_aro3[1], posh_aro3[2], index_atom])
                    atom.append([i, 9, 0.1562, H_aro3])
                    index_atom +=1
                    
                    posc_aro1 = getposc_aro1(posc_aro6, vecmv)
                    N_aro1 = np.asarray([posc_aro1[0], posc_aro1[1], posc_aro1[2], index_atom])
                    atom.append([i, 10, -0.6768, N_aro1])
                    index_atom +=1
                    
                    posc_aro2 = getposc_aro2(posc_aro1, vecd)
                    C_aro2 = np.asarray([posc_aro2[0], posc_aro2[1], posc_aro2[2], index_atom])
                    atom.append([i, 8, 0.4742, C_aro2])
                    index_atom +=1
                    
                    posh_aro2 = getposh_aro2(posc_aro2, vecmv)
                    H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                    atom.append([i, 9, 0.0132, H_aro2])
                    index_atom +=1
                    
                    #Bonding information
                    bonding.append([7, int(C_ali3[3]), int(H_ali4[3])])
                    bonding.append([7, int(C_ali3[3]), int(H_ali5[3])])
                    bonding.append([6, int(C_ali3[3]), int(C_ali4[3])])
                    bonding.append([7, int(C_ali4[3]), int(H_ali3[3])])
                    #bonding.append([1, int(C_ali4[3]), int(C_ali5[3])])
                    #bonding.append([2, int(C_ali5[3]), int(H_ali2[3])])
                    #bonding.append([2, int(C_ali5[3]), int(H_ali1[3])])
                    #bonding.append([2, int(C_ali5[3]), int(H_ali6[3])])
                    bonding.append([7, int(C_ali4[3]), int(H_ali2[3])])
                    bonding.append([11, int(C_ali4[3]), int(C_aro6[3])])
                    bonding.append([8, int(C_aro5[3]), int(C_aro6[3])])
                    bonding.append([10, int(C_aro5[3]), int(H_aro5[3])])
                    bonding.append([8, int(C_aro4[3]), int(C_aro5[3])])
                    bonding.append([10, int(C_aro4[3]), int(H_aro4[3])])
                    bonding.append([8, int(C_aro3[3]), int(C_aro4[3])])
                    bonding.append([10, int(C_aro3[3]), int(H_aro3[3])])
                    bonding.append([9, int(N_aro1[3]), int(C_aro6[3])])
                    bonding.append([9, int(C_aro2[3]), int(N_aro1[3])])
                    bonding.append([8, int(C_aro2[3]), int(C_aro3[3])])
                    bonding.append([10, int(C_aro2[3]), int(H_aro2[3])])
                
            ##      Storing angle information of monomer
                    angle.append([7, int(H_ali4[3]), int(C_ali3[3]), int(H_ali5[3])])
                    angle.append([7, int(H_ali3[3]), int(C_ali4[3]), int(H_ali2[3])])
                    #angle.append([2, int(H_ali2[3]), int(C_ali5[3]), int(H_ali1[3])])
                    #angle.append([2, int(H_ali2[3]), int(C_ali5[3]), int(H_ali6[3])])
                    #angle.append([2, int(H_ali6[3]), int(C_ali5[3]), int(H_ali1[3])])
                    angle.append([6, int(H_ali5[3]), int(C_ali3[3]), int(C_ali4[3])])
                    angle.append([6, int(H_ali3[3]), int(C_ali4[3]), int(C_ali3[3])])
                    angle.append([6, int(H_ali2[3]), int(C_ali4[3]), int(C_ali3[3])])
                    #angle.append([1, int(H_ali3[3]), int(C_ali4[3]), int(C_ali5[3])])
                    angle.append([6, int(H_ali4[3]), int(C_ali3[3]), int(C_ali4[3])])
                    #angle.append([1, int(H_ali2[3]), int(C_ali5[3]), int(C_ali4[3])])
                    #angle.append([1, int(H_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                    #angle.append([1, int(H_ali1[3]), int(C_ali5[3]), int(C_ali4[3])])
                    angle.append([6, int(C_aro6[3]), int(C_ali4[3]), int(H_ali3[3])])
                    angle.append([6, int(C_aro6[3]), int(C_ali4[3]), int(H_ali2[3])])
                    #angle.append([3, int(C_ali3[3]), int(C_ali4[3]), int(C_ali5[3])])
                    angle.append([8, int(C_ali3[3]), int(C_ali4[3]), int(C_aro6[3])])
                    angle.append([9, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([10, int(C_ali4[3]), int(C_aro6[3]), int(N_aro1[3])])
                    angle.append([9, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    angle.append([9, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    angle.append([9, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    angle.append([10, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    angle.append([11, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    angle.append([10, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([12, int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([12, int(C_aro4[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([12, int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([12, int(C_aro3[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([12, int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([12, int(C_aro2[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([13, int(N_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                    angle.append([12, int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    backbone_dihedral.append(int(C_ali3[3]))
                    backbone_dihedral.append(int(C_ali4[3]))
                    backbone_dihedral.append(int(H_ali3[3]))
                    backbone_dihedral.append(int(H_ali2[3]))
                else:
                    #print "func"
                    number_ions +=1
                    
                    x = x + dx
                    y = y + dy
                    z = z + dz
                    
                    posc_ali1 = np.asarray([x,y,z])
                    C_ali3 = np.asarray([x,y,z, index_atom])
                    atom.append([i, 11, 0.0, C_ali3])
                    index_atom +=1
                    
                    C_ali1 = C_ali3
                                
                    posh_ali4 = getposh_ali2(posc_ali1, vecv)
                    H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali4])
                    index_atom +=1
                    
                    posh_ali5 = getposh_ali3(posc_ali1, vecv)
                    H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali5])
                    index_atom +=1
                    
                    #print index_atom,"last monomer"
                    with open("lastmonomer.txt",'a') as outa:
                        outa.write(str(index_atom)+'\n')
                    
                    posc_ali2 = getposc_ali2(posc_ali1, vecm)
                    C_ali4 = np.asarray([posc_ali2[0], posc_ali2[1], posc_ali2[2],index_atom])
                    atom.append([i, 11, 0.0, C_ali4])
                    index_atom +=1
                    
                    C_ali2 = C_ali4
                    
                    posh_ali4 = getposh_ali4(posc_ali2,vecv)
                    H_ali3 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali3])
                    index_atom +=1
                    
                    posh_ali2 = getposh_ali1(posc_ali2,-vec)
                    H_ali2 = np.asarray([posh_ali2[0], posh_ali2[1], posh_ali2[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali2])
                    index_atom +=1
                    
                    posc_aro6 = getposc_aro6(posc_ali2,vecd)
                    C_aro6 = np.asarray([posc_aro6[0], posc_aro6[1], posc_aro6[2], index_atom])
                    atom.append([i, 13, 0.05697, C_aro6])
                    index_atom +=1
                    
                    posc_aro5 = getposc_aro5(posc_aro6,vecmrv)
                    C_aro5 = np.asarray([posc_aro5[0], posc_aro5[1], posc_aro5[2], index_atom])
                    atom.append([i, 13, -0.23303, C_aro5])
                    index_atom +=1
                    
                    posh_aro5 = getposh_aro5(posc_aro5, vecmv)
                    H_aro5 = np.asarray([posh_aro5[0], posh_aro5[1], posh_aro5[2], index_atom])
                    atom.append([i, 14, 0.22957, H_aro5])
                    index_atom +=1
                    
                    posc_aro4 = getposc_aro4(posc_aro5,vecd)
                    C_aro4 = np.asarray([posc_aro4[0], posc_aro4[1], posc_aro4[2], index_atom])
                    atom.append([i, 13, 0.14767, C_aro4])
                    index_atom +=1
                    
                    posh_aro4 = getposh_aro4(posc_aro4,vecmrv)
                    H_aro4 = np.asarray([posh_aro4[0], posh_aro4[1], posh_aro4[2], index_atom])
                    atom.append([i, 14, 0.17877, H_aro4])
                    index_atom +=1
                    
                    posc_aro3 = getposc_aro3(posc_aro6,vecd)
                    C_aro3 = np.asarray([posc_aro3[0], posc_aro3[1], posc_aro3[2], index_atom])
                    atom.append([i, 13, -0.23303, C_aro3])
                    index_atom +=1
                    
                    posh_aro3 = getposh_aro3(posc_aro3,vecd)
                    H_aro3 = np.asarray([posh_aro3[0], posh_aro3[1], posh_aro3[2], index_atom])
                    atom.append([i, 14, 0.22957, H_aro3])
                    index_atom +=1
                    
                    posc_aro1 = getposc_aro1(posc_aro6, vecmv)
                    N_aro1 = np.asarray([posc_aro1[0], posc_aro1[1], posc_aro1[2], index_atom])
                    atom.append([i, 15, 0.19967, N_aro1])
                    index_atom +=1
                    
                    posc_met1 = getposh_aro1(posc_aro1, vecmrv)
                    C_met1 = np.asarray([posc_met1[0], posc_met1[1],posc_met1[2], index_atom])
                    atom.append([i, 16, -0.3979, C_met1])
                    index_atom +=1
                    
                    posh_met1 = getposh_ali1(posc_met1, vecmrv)
                    H_met1 = np.asarray([posh_met1[0], posh_met1[1], posh_met1[2], index_atom])
                    atom.append([i, 17, 0.1916, H_met1])
                    index_atom +=1
                    
                    posh_met2 = getposh_ali2(posc_met1, vecv)
                    H_met2 = np.asarray([posh_met2[0], posh_met2[1], posh_met2[2], index_atom])
                    atom.append([i, 17, 0.1916, H_met2])
                    index_atom +=1
                    
                    posh_met3 = getposh_ali3(posc_met1, vecv)
                    H_met3 = np.asarray([posh_met3[0], posh_met3[1], posh_met3[2], index_atom])
                    atom.append([i, 17, 0.1916, H_met3])
                    index_atom +=1
                    
                    posc_aro2 = getposc_aro2(posc_aro1, vecd)
                    C_aro2 = np.asarray([posc_aro2[0], posc_aro2[1], posc_aro2[2], index_atom])
                    atom.append([i, 13, 0.05697, C_aro2])
                    index_atom +=1
                    
                    posh_aro2 = getposh_aro2(posc_aro2, vecmv)
                    H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                    atom.append([i, 14, 0.18997, H_aro2])
                    index_atom +=1 
                    
                    #Bonding information
                    bonding.append([13, int(C_ali3[3]), int(H_ali4[3])])
                    bonding.append([13, int(C_ali3[3]), int(H_ali5[3])])
                    bonding.append([12, int(C_ali3[3]), int(C_ali4[3])])
                    bonding.append([13, int(C_ali4[3]), int(H_ali3[3])])
                    bonding.append([13, int(C_ali4[3]), int(H_ali2[3])])
                    #bonding.append([1, int(C_ali4[3]), int(C_ali5[3])])
                    #bonding.append([2, int(C_ali5[3]), int(H_ali2[3])])
                    #bonding.append([2, int(C_ali5[3]), int(H_ali1[3])])
                    #bonding.append([2, int(C_ali5[3]), int(H_ali6[3])])
                    
                    bonding.append([17, int(C_ali4[3]), int(C_aro6[3])])
                    bonding.append([14, int(C_aro5[3]), int(C_aro6[3])])
                    bonding.append([16, int(C_aro5[3]), int(H_aro5[3])])
                    bonding.append([14, int(C_aro4[3]), int(C_aro5[3])])
                    bonding.append([16, int(C_aro4[3]), int(H_aro4[3])])
                    bonding.append([14, int(C_aro3[3]), int(C_aro4[3])])
                    bonding.append([16, int(C_aro3[3]), int(H_aro3[3])])
                    bonding.append([15, int(N_aro1[3]), int(C_aro6[3])])
                    bonding.append([15, int(C_aro2[3]), int(N_aro1[3])])
                    bonding.append([14, int(C_aro2[3]), int(C_aro3[3])])
                    bonding.append([16, int(C_aro2[3]), int(H_aro2[3])])
                    ##Due to methylation
                    bonding.append([18, int(N_aro1[3]), int(C_met1[3])])
                    bonding.append([19, int(C_met1[3]), int(H_met1[3])])
                    bonding.append([19, int(C_met1[3]), int(H_met2[3])])
                    bonding.append([19, int(C_met1[3]), int(H_met3[3])])
                    
                    #Storing angle information of monomer
                    angle.append([15, int(H_ali4[3]), int(C_ali3[3]), int(H_ali5[3])])
                    angle.append([15, int(H_ali2[3]), int(C_ali4[3]), int(H_ali3[3])])
                    #angle.append([2, int(H_ali2[3]), int(C_ali5[3]), int(H_ali1[3])])
                    #angle.append([2, int(H_ali2[3]), int(C_ali5[3]), int(H_ali6[3])])
                    #angle.append([2, int(H_ali6[3]), int(C_ali5[3]), int(H_ali1[3])])
                    angle.append([14, int(H_ali5[3]), int(C_ali3[3]), int(C_ali4[3])])
                    angle.append([14, int(H_ali3[3]), int(C_ali4[3]), int(C_ali3[3])])
                    angle.append([14, int(H_ali2[3]), int(C_ali4[3]), int(C_ali3[3])])
                    #angle.append([1, int(H_ali3[3]), int(C_ali4[3]), int(C_ali5[3])])
                    angle.append([14, int(H_ali4[3]), int(C_ali3[3]), int(C_ali4[3])])
                    #angle.append([1, int(H_ali2[3]), int(C_ali5[3]), int(C_ali4[3])])
                    #angle.append([1, int(H_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                    #angle.append([1, int(H_ali1[3]), int(C_ali5[3]), int(C_ali4[3])])
                    angle.append([14, int(C_aro6[3]), int(C_ali4[3]), int(H_ali3[3])])
                    angle.append([14, int(C_aro6[3]), int(C_ali4[3]), int(H_ali2[3])])
                    #angle.append([3, int(C_ali3[3]), int(C_ali4[3]), int(C_ali5[3])])
                    angle.append([16, int(C_ali3[3]), int(C_ali4[3]), int(C_aro6[3])])
                    angle.append([17, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([18, int(C_ali4[3]), int(C_aro6[3]), int(N_aro1[3])])
                    angle.append([17, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    angle.append([17, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    angle.append([17, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    angle.append([18, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    angle.append([19, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    angle.append([18, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([20, int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([20, int(C_aro4[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([20, int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([20, int(C_aro3[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([20, int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([20, int(C_aro2[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([21, int(N_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                    angle.append([20, int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    ##Due to methylation
                    angle.append([19, int(C_aro6[3]), int(N_aro1[3]), int(C_met1[3])])
                    angle.append([19, int(C_aro2[3]), int(N_aro1[3]), int(C_met1[3])])
                    angle.append([22, int(N_aro1[3]), int(C_met1[3]), int(H_met1[3])])
                    angle.append([22, int(N_aro1[3]), int(C_met1[3]), int(H_met2[3])])
                    angle.append([22, int(N_aro1[3]), int(C_met1[3]), int(H_met3[3])])
                    angle.append([23, int(H_met1[3]), int(C_met1[3]), int(H_met2[3])])
                    angle.append([23, int(H_met1[3]), int(C_met1[3]), int(H_met3[3])])
                    angle.append([23, int(H_met2[3]), int(C_met1[3]), int(H_met3[3])])
                    
                    #Storing dihedrals
                    #dihedral.append([1, int(H_ali2[3]), int(C_ali5[3]), int(C_ali4[3]), int(C_ali3[3])])
                    #dihedral.append([1, int(H_ali1[3]), int(C_ali5[3]), int(C_ali4[3]), int(C_ali3[3])])
                    #dihedral.append([1, int(H_ali6[3]), int(C_ali5[3]), int(C_ali4[3]), int(C_ali3[3])])
                    backbone_dihedral.append(int(C_ali3[3]))
                    backbone_dihedral.append(int(C_ali4[3]))
                    backbone_dihedral.append(int(H_ali2[3]))
                    backbone_dihedral.append(int(H_ali3[3]))
                    dihedral.append([2, int(C_ali2[3]), int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    dihedral.append([2, int(C_ali2[3]), int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    dihedral.append([2, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    dihedral.append([2, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    dihedral.append([2, int(H_aro5[3]), int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    dihedral.append([2, int(H_aro5[3]), int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    dihedral.append([2, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    dihedral.append([2, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    dihedral.append([2, int(H_aro4[3]), int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    dihedral.append([2, int(H_aro4[3]), int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    dihedral.append([2, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    dihedral.append([2, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    dihedral.append([2, int(H_aro3[3]), int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    dihedral.append([2, int(H_aro3[3]), int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    dihedral.append([3, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([4, int(H_aro2[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    dihedral.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_ali2[3])])
                    dihedral.append([2, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    dihedral.append([2, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    ##Due to methylation
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_ali2[3])])
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3]), int(C_aro3[3])])
                    dihedral.append([5, int(H_met1[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([5, int(H_met2[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([5, int(H_met3[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([5, int(H_met1[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3])])
                    dihedral.append([5, int(H_met2[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3])])
                    dihedral.append([5, int(H_met3[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3])])
                    
                    ## Improper information
                    improper.append([2, int(C_aro6[3]), int(N_aro1[3]), int(C_aro5[3]), int(C_ali2[3])])
                    improper.append([2, int(C_aro5[3]), int(C_aro6[3]), int(C_aro4[3]), int(H_aro5[3])])
                    improper.append([2, int(C_aro4[3]), int(C_aro5[3]), int(C_aro3[3]), int(H_aro4[3])])
                    improper.append([2, int(C_aro3[3]), int(C_aro4[3]), int(C_aro2[3]), int(H_aro3[3])])
                    improper.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro3[3]), int(H_aro2[3])])
                    improper.append([2, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    improper.append([2, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    improper.append([2, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    improper.append([2, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    improper.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    improper.append([2, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    #Due to methylation
                    improper.append([2, int(C_aro6[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_met1[3])])
                    
                    #Adding backbone dihedrals
                    backbone_dihedral = np.asarray(backbone_dihedral)
                    for b in range(len(backbone_dihedral)-4):
                    #dihedral.append([1, backbone_dihedral[i], backbone_dihedral[i+1], backbone_dihedral[i+2], backbone_dihedral[i+3]])
                            if b ==len(backbone_dihedral)-4:
                                dihedral.append([1, backbone_dihedral[b], backbone_dihedral[b+1], backbone_dihedral[b+2], backbone_dihedral[b+4]])
                                dihedral.append([1, backbone_dihedral[b], backbone_dihedral[b+1], backbone_dihedral[b+2], backbone_dihedral[b+3]])
                            
                            else:
                                dihedral.append([1, backbone_dihedral[b], backbone_dihedral[b+1], backbone_dihedral[b+2], backbone_dihedral[b+3]])   
            else:
                if(j==0):
                    #print "for internal monomer,j=0", i*nB+k+1
                   # print " for internal monomer, j = 0", j ,k,i
                    if(k == nB-1):
                        vec = np.asarray([t1[i*nB+k][0]-t1[i*nB+k-1][0],t1[i*nB+k][1]-t1[i*nB+k-1][1],t1[i*nB+k][2]-t1[i*nB+k-1][2]])
                    else:
                        vec = np.asarray([t1[i*nB+k+1][0]-t1[i*nB+k][0],t1[i*nB+k+1][1]-t1[i*nB+k][1],t1[i*nB+k+1][2]-t1[i*nB+k][2]])
                    lentv = lent(vec)
                    #the vector which is the direction of two bead
                    vec = vec/lent(vec)
                    #a unit vector which is vertical t vec
                   # print vec
                    vecv = dirv(vec)
                    vecd = dirvd(vec, vecv)
                    vecm = 1.732*vec + vecd
                    vecm = vecm/lent(vecm)
                    vecmr = vecd - 1.732*vec
                    vecmr = vecmr/lent(vecmr)
                    vecmv = 1.732*vecv +vecd
                    vecmv = vecmv/lent(vecmv)
                    vecmrv = vecd - 1.732*vecv
                    vecmrv = vecmrv/lent(vecmrv)
                    dx = vec[0]*lentv/mpb
                    dy = vec[1]*lentv/mpb
                    dz = vec[2]*lentv/mpb 
                    x = t1[i*nA+k][0]
                    y = t1[i*nA+k][1]
                    z = t1[i*nA+k][2] 
                else:
                   # print " for internal monomer", j ,k, i*nB+k
                    x = x + dx
                    y = y + dy
                    z = z + dz
                funct = random.random()
                #print funct 
                #print func, "func"   
                if (funct > func):
                    #print " no func"
                    posc_ali1 = np.asarray([x,y,z])
                    C_ali3 = np.asarray([x,y,z, index_atom])
                    atom.append([i, 6, 0.0, C_ali3])
                    index_atom +=1
                    
                    
                    C_ali1 = C_ali3
                                
                    posh_ali4 = getposh_ali2(posc_ali1, vecv)
                    H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali4])
                    index_atom +=1
                    
                    posh_ali5 = getposh_ali3(posc_ali1, vecv)
                    H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali5])
                    index_atom +=1
                    
                    posc_ali2 = getposc_ali2(posc_ali1, vecm)
                    C_ali4 = np.asarray([posc_ali2[0], posc_ali2[1], posc_ali2[2],index_atom])
                    atom.append([i, 6, 0.0, C_ali4])
                    index_atom +=1
                    
                    C_ali2 = C_ali4
                    
                    posh_ali4 = getposh_ali4(posc_ali2,vecv)
                    H_ali3 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 7, 0.0, H_ali3])
                    index_atom +=1
                    
                    posc_aro6 = getposc_aro6(posc_ali2,vecd)
                    C_aro6 = np.asarray([posc_aro6[0], posc_aro6[1], posc_aro6[2], index_atom])
                    atom.append([i, 8, 0.4742, C_aro6])
                    index_atom +=1
                    
                    posc_aro5 = getposc_aro5(posc_aro6,vecmrv)
                    C_aro5 = np.asarray([posc_aro5[0], posc_aro5[1], posc_aro5[2], index_atom])
                    atom.append([i, 8, -0.4458, C_aro5])
                    index_atom +=1
                    
                    posh_aro5 = getposh_aro5(posc_aro5, vecmv)
                    H_aro5 = np.asarray([posh_aro5[0], posh_aro5[1], posh_aro5[2], index_atom])
                    atom.append([i, 9, 0.1562, H_aro5])
                    index_atom +=1
                    
                    posc_aro4 = getposc_aro4(posc_aro5,vecd)
                    C_aro4 = np.asarray([posc_aro4[0], posc_aro4[1], posc_aro4[2], index_atom])
                    atom.append([i, 8, 0.2282, C_aro4])
                    index_atom +=1
                    
                    posh_aro4 = getposh_aro4(posc_aro4,vecmrv)
                    H_aro4 = np.asarray([posh_aro4[0], posh_aro4[1], posh_aro4[2], index_atom])
                    atom.append([i, 9, 0.0662, H_aro4])
                    index_atom +=1
                    
                    posc_aro3 = getposc_aro3(posc_aro6,vecd)
                    C_aro3 = np.asarray([posc_aro3[0], posc_aro3[1], posc_aro3[2], index_atom])
                    atom.append([i, 8, -0.4458, C_aro3])
                    index_atom +=1
                    
                    posh_aro3 = getposh_aro3(posc_aro3,vecd)
                    H_aro3 = np.asarray([posh_aro3[0], posh_aro3[1], posh_aro3[2], index_atom])
                    atom.append([i, 9, 0.1562, H_aro3])
                    index_atom +=1
                    
                    posc_aro1 = getposc_aro1(posc_aro6, vecmv)
                    N_aro1 = np.asarray([posc_aro1[0], posc_aro1[1], posc_aro1[2], index_atom])
                    atom.append([i, 10, -0.6768, N_aro1])
                    index_atom +=1
                    
                    posc_aro2 = getposc_aro2(posc_aro1, vecd)
                    C_aro2 = np.asarray([posc_aro2[0], posc_aro2[1], posc_aro2[2], index_atom])
                    atom.append([i, 8, 0.4742, C_aro2])
                    index_atom +=1
                    
                    posh_aro2 = getposh_aro2(posc_aro2, vecmv)
                    H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                    atom.append([i, 9, 0.0132, H_aro2])
                    index_atom +=1
                    
                    ### Bonding information for monomer
                    bonding.append([7, int(C_ali3[3]), int(H_ali4[3])])
                    bonding.append([7, int(C_ali3[3]), int(H_ali5[3])])
                    bonding.append([6, int(C_ali3[3]), int(C_ali4[3])])
                    bonding.append([7, int(C_ali4[3]), int(H_ali3[3])])
                    bonding.append([11, int(C_ali4[3]), int(C_aro6[3])])
                    bonding.append([8, int(C_aro5[3]), int(C_aro6[3])])
                    bonding.append([10, int(C_aro5[3]), int(H_aro5[3])])
                    bonding.append([8, int(C_aro4[3]), int(C_aro5[3])])
                    bonding.append([10, int(C_aro4[3]), int(H_aro4[3])])
                    bonding.append([8, int(C_aro3[3]), int(C_aro4[3])])
                    bonding.append([10, int(C_aro3[3]), int(H_aro3[3])])
                    bonding.append([9, int(N_aro1[3]), int(C_aro6[3])])
                    bonding.append([9, int(C_aro2[3]), int(N_aro1[3])])
                    bonding.append([8, int(C_aro2[3]), int(C_aro3[3])])
                    bonding.append([10, int(C_aro2[3]), int(H_aro2[3])])
                    bonding.append([6, int(C_ali4[3]), int(index_atom)]) #bonding to next monomer
                 
                    
                    
                    
        ##Storing angle information of monomer
                    angle.append([7, int(H_ali4[3]), int(C_ali3[3]), int(H_ali5[3])])
                    angle.append([6, int(H_ali5[3]), int(C_ali3[3]), int(C_ali4[3])])
                    angle.append([6, int(H_ali3[3]), int(C_ali4[3]), int(C_ali3[3])])
                    angle.append([6, int(H_ali4[3]), int(C_ali3[3]), int(C_ali4[3])])
                    angle.append([6, int(H_ali3[3]), int(C_ali4[3]), int(index_atom)])
                    angle.append([6, int(index_atom+1), int(index_atom), int(C_ali4[3])])
                    angle.append([6, int(index_atom+2), int(index_atom), int(C_ali4[3])])
                    angle.append([6, int(C_aro6[3]), int(C_ali4[3]), int(H_ali3[3])])
                    angle.append([8, int(C_ali3[3]), int(C_ali4[3]), int(index_atom)])
                    angle.append([8, int(C_ali4[3]), int(index_atom), int(index_atom+3)])
                    angle.append([8, int(C_ali3[3]), int(C_ali4[3]), int(C_aro6[3])])
                    angle.append([8, int(index_atom), int(C_ali4[3]), int(C_aro6[3])])
                    angle.append([9, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([10, int(C_ali4[3]), int(C_aro6[3]), int(N_aro1[3])])
                    angle.append([9, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    angle.append([9, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    angle.append([9, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    angle.append([10, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    angle.append([11, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    angle.append([10, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([12, int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([12, int(C_aro4[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([12, int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([12, int(C_aro3[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([12, int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([12, int(C_aro2[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([13, int(N_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                    angle.append([12, int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    
            ##Storing dihedral information of monomer
                    #dihedral.append([1, int(C_ali3[3]), int(C_ali4[3]), int(index_atom), int(index_atom+3)])
                    
                    #if j+1==num_int_monomers_P2VP: ##penultimate monomer has different dihedral
                    #   dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+5)])
                    #   dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+4)])
                    #else:  
                        #dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+15)])
                    backbone_dihedral.append(int(C_ali3[3]))
                    backbone_dihedral.append(int(C_ali4[3]))
                else:
                    
                                        
                    number_ions +=1
                    
                    posc_ali1 = np.asarray([x,y,z])
                    C_ali3 = np.asarray([x,y,z, index_atom])
                    atom.append([i, 11, 0.0, C_ali3])
                    index_atom +=1
                    
                    C_ali1 = C_ali3
                                
                    posh_ali4 = getposh_ali2(posc_ali1, vecv)
                    H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali4])
                    index_atom +=1
                    
                    posh_ali5 = getposh_ali3(posc_ali1, vecv)
                    H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali5])
                    index_atom +=1
                    
                    posc_ali2 = getposc_ali2(posc_ali1, vecm)
                    C_ali4 = np.asarray([posc_ali2[0], posc_ali2[1], posc_ali2[2],index_atom])
                    atom.append([i, 11, 0.0, C_ali4])
                    index_atom +=1
                    
                    C_ali2 = C_ali4
                    
                    posh_ali4 = getposh_ali4(posc_ali2,vecv)
                    H_ali3 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                    atom.append([i, 12, 0.0, H_ali3])
                    index_atom +=1
                   
                    posc_aro6 = getposc_aro6(posc_ali2,vecd)
                    C_aro6 = np.asarray([posc_aro6[0], posc_aro6[1], posc_aro6[2], index_atom])
                    atom.append([i, 13, 0.05697, C_aro6])
                    index_atom +=1
                    
                    posc_aro5 = getposc_aro5(posc_aro6,vecmrv)
                    C_aro5 = np.asarray([posc_aro5[0], posc_aro5[1], posc_aro5[2], index_atom])
                    atom.append([i, 13, -0.23303, C_aro5])
                    index_atom +=1
                    
                    posh_aro5 = getposh_aro5(posc_aro5, vecmv)
                    H_aro5 = np.asarray([posh_aro5[0], posh_aro5[1], posh_aro5[2], index_atom])
                    atom.append([i, 14, 0.22957, H_aro5])
                    index_atom +=1
                    
                    posc_aro4 = getposc_aro4(posc_aro5,vecd)
                    C_aro4 = np.asarray([posc_aro4[0], posc_aro4[1], posc_aro4[2], index_atom])
                    atom.append([i, 13, 0.14767, C_aro4])
                    index_atom +=1
                    
                    posh_aro4 = getposh_aro4(posc_aro4,vecmrv)
                    H_aro4 = np.asarray([posh_aro4[0], posh_aro4[1], posh_aro4[2], index_atom])
                    atom.append([i, 14, 0.17877, H_aro4])
                    index_atom +=1
                    
                    posc_aro3 = getposc_aro3(posc_aro6,vecd)
                    C_aro3 = np.asarray([posc_aro3[0], posc_aro3[1], posc_aro3[2], index_atom])
                    atom.append([i, 13, -0.23303, C_aro3])
                    index_atom +=1
                    
                    posh_aro3 = getposh_aro3(posc_aro3,vecd)
                    H_aro3 = np.asarray([posh_aro3[0], posh_aro3[1], posh_aro3[2], index_atom])
                    atom.append([i, 14, 0.22957, H_aro3])
                    index_atom +=1
                    
                    posc_aro1 = getposc_aro1(posc_aro6, vecmv)
                    N_aro1 = np.asarray([posc_aro1[0], posc_aro1[1], posc_aro1[2], index_atom])
                    atom.append([i, 15, 0.19967, N_aro1])
                    index_atom +=1
                    
                    posc_met1 = getposh_aro1(posc_aro1, vecmrv)
                    C_met1 = np.asarray([posc_met1[0], posc_met1[1],posc_met1[2], index_atom])
                    atom.append([i, 16, -0.3979, C_met1])
                    index_atom +=1
                    
                    posh_met1 = getposh_ali1(posc_met1, vecmrv)
                    H_met1 = np.asarray([posh_met1[0], posh_met1[1], posh_met1[2], index_atom])
                    atom.append([i, 17, 0.1916, H_met1])
                    index_atom +=1
                    
                    posh_met2 = getposh_ali2(posc_met1, vecv)
                    H_met2 = np.asarray([posh_met2[0], posh_met2[1], posh_met2[2], index_atom])
                    atom.append([i, 17, 0.1916, H_met2])
                    index_atom +=1
                    
                    posh_met3 = getposh_ali3(posc_met1, vecv)
                    H_met3 = np.asarray([posh_met3[0], posh_met3[1], posh_met3[2], index_atom])
                    atom.append([i, 17, 0.1916, H_met3])
                    index_atom +=1
                    
                    posc_aro2 = getposc_aro2(posc_aro1, vecd)
                    C_aro2 = np.asarray([posc_aro2[0], posc_aro2[1], posc_aro2[2], index_atom])
                    atom.append([i, 13, 0.05697, C_aro2])
                    index_atom +=1
                    
                    posh_aro2 = getposh_aro2(posc_aro2, vecmv)
                    H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                    atom.append([i, 14, 0.18997, H_aro2])
                    index_atom +=1 
                    
                    #Storing information for monomer
                    bonding.append([13, int(C_ali3[3]), int(H_ali4[3])])
                    bonding.append([13, int(C_ali3[3]), int(H_ali5[3])])
                    bonding.append([12, int(C_ali3[3]), int(C_ali4[3])])
                    bonding.append([13, int(C_ali4[3]), int(H_ali3[3])])
                    bonding.append([17, int(C_ali4[3]), int(C_aro6[3])])
                    bonding.append([14, int(C_aro5[3]), int(C_aro6[3])])
                    bonding.append([16, int(C_aro5[3]), int(H_aro5[3])])
                    bonding.append([14, int(C_aro4[3]), int(C_aro5[3])])
                    bonding.append([16, int(C_aro4[3]), int(H_aro4[3])])
                    bonding.append([14, int(C_aro3[3]), int(C_aro4[3])])
                    bonding.append([16, int(C_aro3[3]), int(H_aro3[3])])
                    bonding.append([15, int(N_aro1[3]), int(C_aro6[3])])
                    bonding.append([15, int(C_aro2[3]), int(N_aro1[3])])
                    bonding.append([14, int(C_aro2[3]), int(C_aro3[3])])
                    bonding.append([16, int(C_aro2[3]), int(H_aro2[3])])
                    bonding.append([12, int(C_ali2[3]), int(index_atom)]) #bonding to next monomer
                    ##Due to methylation
                    bonding.append([18, int(N_aro1[3]), int(C_met1[3])])
                    bonding.append([19, int(C_met1[3]), int(H_met1[3])])
                    bonding.append([19, int(C_met1[3]), int(H_met2[3])])
                    bonding.append([19, int(C_met1[3]), int(H_met3[3])])
                    
                    #Storing angle information of monomer
                    angle.append([15, int(H_ali4[3]), int(C_ali3[3]), int(H_ali5[3])])
                    angle.append([14, int(H_ali5[3]), int(C_ali3[3]), int(C_ali4[3])])
                    angle.append([14, int(H_ali3[3]), int(C_ali4[3]), int(C_ali3[3])])
                    angle.append([14, int(H_ali4[3]), int(C_ali3[3]), int(C_ali4[3])])
                    angle.append([14, int(H_ali3[3]), int(C_ali4[3]), int(index_atom)])
                    angle.append([14, int(index_atom+1), int(index_atom), int(C_ali4[3])])
                    angle.append([14, int(index_atom+2), int(index_atom), int(C_ali4[3])])
                    angle.append([14, int(C_aro6[3]), int(C_ali4[3]), int(H_ali3[3])])
                    angle.append([16, int(C_ali3[3]), int(C_ali4[3]), int(index_atom)])
                    angle.append([16, int(C_ali4[3]), int(index_atom), int(index_atom+3)])
                    angle.append([16, int(C_ali3[3]), int(C_ali4[3]), int(C_aro6[3])])
                    angle.append([16, int(index_atom), int(C_ali4[3]), int(C_aro6[3])])
                    angle.append([17, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([18, int(C_ali4[3]), int(C_aro6[3]), int(N_aro1[3])])
                    angle.append([17, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    angle.append([17, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    angle.append([17, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    angle.append([18, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    angle.append([19, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    angle.append([18, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    angle.append([20, int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([20, int(C_aro4[3]), int(C_aro5[3]), int(H_aro5[3])])
                    angle.append([20, int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([20, int(C_aro3[3]), int(C_aro4[3]), int(H_aro4[3])])
                    angle.append([20, int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([20, int(C_aro2[3]), int(C_aro3[3]), int(H_aro3[3])])
                    angle.append([21, int(N_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                    angle.append([20, int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    ##Due to methylation
                    angle.append([19, int(C_aro6[3]), int(N_aro1[3]), int(C_met1[3])])
                    angle.append([19, int(C_aro2[3]), int(N_aro1[3]), int(C_met1[3])])
                    angle.append([22, int(N_aro1[3]), int(C_met1[3]), int(H_met1[3])])
                    angle.append([22, int(N_aro1[3]), int(C_met1[3]), int(H_met2[3])])
                    angle.append([22, int(N_aro1[3]), int(C_met1[3]), int(H_met3[3])])
                    angle.append([23, int(H_met1[3]), int(C_met1[3]), int(H_met2[3])])
                    angle.append([23, int(H_met1[3]), int(C_met1[3]), int(H_met3[3])])
                    angle.append([23, int(H_met2[3]), int(C_met1[3]), int(H_met3[3])])
                    
                    #Storing dihedral information of monomer
                    #dihedral.append([1, int(C_ali3[3]), int(C_ali4[3]), int(index_atom), int(index_atom+3)])
                    #
                    #if j+1==num_int_monomers: ##penultimate monomer has different dihedral
                    #   dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+5)])
                    #   dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+4)])
                    #else:  
                    #   dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+19)])
                    backbone_dihedral.append(int(C_ali3[3]))
                    backbone_dihedral.append(int(C_ali4[3]))
                    
                    dihedral.append([2, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    dihedral.append([2, int(C_ali4[3]), int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    dihedral.append([2, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    dihedral.append([2, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    dihedral.append([2, int(H_aro5[3]), int(C_aro5[3]), int(C_aro4[3]), int(H_aro4[3])])
                    dihedral.append([2, int(H_aro5[3]), int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    dihedral.append([2, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    dihedral.append([2, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    dihedral.append([2, int(H_aro4[3]), int(C_aro4[3]), int(C_aro3[3]), int(H_aro3[3])])
                    dihedral.append([2, int(H_aro4[3]), int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    dihedral.append([2, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    dihedral.append([2, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    dihedral.append([2, int(H_aro3[3]), int(C_aro3[3]), int(C_aro2[3]), int(H_aro2[3])])
                    dihedral.append([2, int(H_aro3[3]), int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    dihedral.append([3, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([4, int(H_aro2[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    dihedral.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_ali2[3])])
                    dihedral.append([2, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3]), int(H_aro5[3])])
                    dihedral.append([2, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    ##Due to methylation
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_ali2[3])])
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3]), int(H_aro2[3])])
                    dihedral.append([4, int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3]), int(C_aro3[3])])
                    dihedral.append([5, int(H_met1[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([5, int(H_met2[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([5, int(H_met3[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro6[3])])
                    dihedral.append([5, int(H_met1[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3])])
                    dihedral.append([5, int(H_met2[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3])])
                    dihedral.append([5, int(H_met3[3]), int(C_met1[3]), int(N_aro1[3]), int(C_aro2[3])])
                    
                    ## Improper information
                    improper.append([2, int(C_aro6[3]), int(N_aro1[3]), int(C_aro5[3]), int(C_ali2[3])])
                    improper.append([2, int(C_aro5[3]), int(C_aro6[3]), int(C_aro4[3]), int(H_aro5[3])])
                    improper.append([2, int(C_aro4[3]), int(C_aro5[3]), int(C_aro3[3]), int(H_aro4[3])])
                    improper.append([2, int(C_aro3[3]), int(C_aro4[3]), int(C_aro2[3]), int(H_aro3[3])])
                    improper.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro3[3]), int(H_aro2[3])])
                    improper.append([2, int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3])])
                    improper.append([2, int(C_aro5[3]), int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3])])
                    improper.append([2, int(C_aro4[3]), int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3])])
                    improper.append([2, int(C_aro3[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3])])
                    improper.append([2, int(C_aro2[3]), int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3])])
                    improper.append([2, int(N_aro1[3]), int(C_aro6[3]), int(C_aro5[3]), int(C_aro4[3])])
                    #Due to methylation
                    improper.append([2, int(C_aro6[3]), int(C_aro2[3]), int(N_aro1[3]), int(C_met1[3])])
                    


total_atoms = len(atom)
total_bonds = len(bonding)
total_angles = len(angle)
total_dihedrals = len(dihedral)
total_impropers = len(improper)


f = open("PS_P2VP_random.data","w")
###Headers
with open("PS_P2VP_random.data", "a") as out_file:
    out_file.write('PS with ' + str(num_chains) +' chains and ' + str(mpb*nA-1) + ' internal PS monomers per chain and ' + str(mpb*nB-1) +' int P2VP monomers per chain'+ 'and ' + str(number_ions) + ' ions' +'\n'
    + '\n' +
    str(total_atoms) + ' ' + 'atoms' + '\n'
    + str(total_bonds) + ' ' + 'bonds' + '\n'
    + str(total_angles) + ' ' + 'angles' + '\n'
    + str(total_dihedrals) + ' ' + 'dihedrals' + '\n'
    + str(total_impropers) + ' ' + 'impropers' + '\n'
    + '\n'
    + """20 atom types
20 bond types
24 angle types
5 dihedral types
2 improper types

0 200 xlo xhi
0 200 ylo yhi
0 200 zlo zhi

Masses

1 12
2 1.008
3 12
4 12
5 1.008
6 12
7 1.008
8 12
9 1.008
10 14.007
11 12
12 1.008
13 12
14 1.008
15 14.007
16 12
17 1.008 
18 126.9
19 15.9994
20 1.008

Bond Coeffs

1 310.44 1.53 #134 1.53
2 317.45 1.51 #158.5 1.51
3 469.47 1.39 #234.5 1.39
4 239.00 1.1  #170 1.1
5 239.00 1.08 #183.5 1.08
6 310.44 1.53 #134 1.53
7 239.00 1.09 #170 1.1
8 469.0 1.40
9 483.0 1.339
10 367.0 1.08
11 317.45 1.51 #158.5 1.51
12 310.44 1.53 #134 1.53
13 239.00 1.09 #170 1.1
14 469.0 1.40
15 483.0 1.339
16 367.0 1.08
17 317.45 1.51 #158.5 1.51
18 337.0 1.481
19 340.0 1.079 
20 1000.0 1.0

Angle Coeffs

1 43.8456 109.45
2 36.6157 109.45
3 57.6362 109.45
4 45.0048 120
5 50.0478 120
6 43.8456 109.45
7 36.6157 109.45
8 57.6362 109.45
9 63 120
10 70 124
11 70 117
12 35 120
13 35 116
14 43.8456 109.45
15 36.6157 109.45
16 57.6362 109.45
17 63 120
18 70 124
19 70 117
20 35 120
21 35 116
22 37.5 109
23 33 109.9 
24 100.0 109.47

Dihedral Coeffs

1 1.43403 -1 3
2 3.625 -1 2 
3 0.15 -1 2 
4 1.5 -1 2
5 1.75 -1 2 

Improper Coeffs

1 harmonic 20.0048 0
2 cvff 1.1 -1 2""" + '\n' + '\n' + 'Atoms' +'\n' +'\n')

## add Atoms

for i in range(len(atom)):
    write_data_atom(i+1, atom[i])
#with open("PS_P2VP_sac.data", "r") as in_file:
#   with open("PS_P2VP_random.data", "a") as out_file:
#       out_file.write('\n' + "Atoms" + '\n' + '\n')
#       all_coords = in_file.readlines()
#       out_file.write(all_coords)


## Adding Bonding information
bonding=np.asarray(bonding) 
with open("PS_P2VP_random.data", "a") as out_file:
    out_file.write('\n'+"Bonds" + '\n' + '\n')
for i in range(len(bonding)):
    write_data_bond(i+1, bonding[i])
## Adding Angle information
angle = np.asarray(angle)
with open("PS_P2VP_random.data", "a") as out_file:
    out_file.write('\n' + "Angles" + '\n' + '\n')
for i in range(len(angle)):
    write_data_angle(i+1, angle[i])
##Adding dihedral information
dihedral = np.asarray(dihedral)
with open("PS_P2VP_random.data","a") as out_file:
    out_file.write('\n' + "Dihedrals" + '\n' + '\n')
for i in range(len(dihedral)):
    write_data_dihedral(i+1, dihedral[i])
##Adding improper information
improper = np.asarray(improper)
with open("PS_P2VP_random.data", "a") as out_file:
    out_file.write('\n' + "Impropers" + '\n' + '\n')
for i in range(len(improper)):
    write_data_improper(i+1, improper[i])
