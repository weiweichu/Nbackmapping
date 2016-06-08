#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  backmappingN2.py
#  
#  Copyright 2016 weiwei <weiwei@xps8700>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

import sys
import math
import pylab
import matplotlib.pyplot as plt
import numpy as np
import random

filename = "config.out"
nA = 8
nB = 8
nAB = nA+nB
density = 0.7
mpb = 8 #monomer per bead
bondlength = 1.5
func = 0
xl = 0
xh = 6
yl = 0
yh = 6
zl = 0
zh = 6

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
    g = np.array([0,1,0])
    g[2] = -a[1]/a[2]
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
    with open("PS_PMMA_50.data", "a") as out_file:
        #out_file.write(str(index_atom) + ' ' + str(atom).strip('[]').strip(',') + '\n')
        out_file.write(str(index_atom)+ ' ' + str(index_chain)+ ' ' + str(type_atom) + ' ' + str(charge) + ' ' + 
        str("{0:.3f}".format(x_coord)).strip('[]') + ' ' + str("{0:.3f}".format(y_coord)).strip('[]') + ' ' +
        str("{0:.3f}".format(z_coord)).strip('[]') + '\n')#+ ' ' + name + ' ' + monomer)

def write_data_bond(index_bond, type_bond):#, atom1, atom2):
    with open("PS_PMMA_50.data", "a") as out_file:
        out_file.write(str(index_bond) + ' ' + str(type_bond).strip('[]').strip(',')+'\n')# + ' ' + str(atom1) + ' ' + str(atom2) + '\n')

def write_data_angle(index_angle, type_angle):#, atom1, atom2, atom3):
    with open("PS_PMMA_50.data", "a") as out_file:
        out_file.write(str(index_angle) + ' ' + str(type_angle).strip('[]').strip(',') + '\n')# + ' ' + str(atom1) + ' ' + str(atom2) + ' ' + str(atom3) + '\n')

def write_data_dihedral(index_dihedral, type_dihedral):#, atom1, atom2, atom3, atom4):
    with open("PS_PMMA_50.data", "a") as out_file:
        out_file.write(str(index_dihedral) + ' ' + str(type_dihedral).strip('[]').strip(',') + '\n')# + ' ' + str(atom1) + ' ' + str(atom2) + ' ' + str(atom3) + ' ' + str(atom4) + '\n')

def write_data_improper(index_improper, type_improper):#, atom1, atom2, atom3, atom4):
    with open("PS_PMMA_50.data", "a") as out_file:
        out_file.write(str(index_improper) + ' ' + str(type_improper).strip('[]').strip(',') + '\n')# + ' ' + str(atom1) + ' ' + str(atom2) + ' ' + str(atom3) + ' ' + str(atom4) + '\n')   
                            
g = open(filename, 'r')
linesa = g.readlines()
perchain = [0,0,0,0]
k=0
chain = 0
d = open('sconfig.xyz', 'w')
i = 2

#select chains according to defined density and in a certain region(xl, xh,
# yl, yh, zl, zh) and rescale 

while i < density*len(linesa):
    #print i
    linesa[i] = linesa[i].split()
    linesa[i] = np.array(map(float, linesa[i]))
    perchain = np.vstack((perchain, linesa[i]))
    k = k+1
    if(k%nAB != 0):
        if(perchain[k][1] > xl and perchain[k][1] < xh and perchain[k][2] > yl and perchain[k][2] < yh and perchain[k][3] > zl and perchain[k][3] < zh):
            i += 1
        else:
            i = i+nAB+1-k
            k = 0
            perchain = [0,0,0,0]
    else:
        d = open('sconfig.xyz','a')
        perchain = perchain[1:]
        for m in range(nAB):
            d.write('{0:8d}{1:15.5f}{2:15.5f}{3:15.5f}\n'.format(int(perchain[m][0]), perchain[m][1], perchain[m][2],perchain[m][3]))
        perchain = [0,0,0,0]
        k = 0
        chain += 1
        i += 1
nPolymer = chain
print 'nPolymer', nPolymer  
g.close()
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
print bondlength, mpb, bond
sf = bondlength * 1.732 * mpb/bond
print "the scale factor is ", sf
bond = 0
k=0
for i in range(len(t0)):
    t0[i] = [t0[i][0]*sf, t0[i][1]*sf, t0[i][2]*sf]
for i in range(len(t1)):
    t1[i] = [t1[i][0]*sf, t1[i][1]*sf, t1[i][2]*sf]
    
for i in range(1,len(t0)):
    if(i%nA != 0):
        bond += lenp(t0[i-1],t0[i])
        k += 1
bond /= k
print "bondlength for scaled coarse grain bead  ", bond
bondm = 1.732*bondlength



#back mapping N polymer chains with nA A beads and nB B beads, one bead 
#will be back mapped to nM monomers

num_chains = nPolymer
num_int_monomers_PS = nA -1
num_int_monomers_PMMA = nB - 1
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
  

for i in range(num_chains):
    #print "chain",i
    backbone_dihedral = []
    for k in range(nA):
        if(k == nA-1):
            mpb = int(lenp(t1[i*nB],t0[i*nA+k])/bondm)
        else:
            mpb = int(lenp(t0[i*nA+k+1],t0[i*nA+k])/bondm)
        for j in range(mpb):
            if(k == 0 and j == 0):
               # print " k==0,j==0", k,j
                #Initial monomer of PS
                vec = np.asarray([t0[i*nA+k+1][0]-t0[i*nA+k][0],t0[i*nA+k+1][1]-t0[i*nA+k][1],t0[i*nA+k+1][2]-t0[i*nA+k][2]])
                #print k, j
               # print "vec", vec
                lentv = lent(vec)
                #the vector which is the direction of two bead
                vec = vec/lent(vec)
                #a unit vector which is vertical t vec
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
                
                posh_aro2 = getposh_aro2(posc_aro2, vecmv)
                H_aro2 = np.asarray([posh_aro2[0], posh_aro2[1], posh_aro2[2], index_atom])
                atom.append([i, 5, 0.115, H_aro2])
                index_atom +=1
                
                ##Storing bonding information of monomer
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
                
            ##Storing dihedral information of monomer
                #dihedral.append([1, int(C_ali3[3]), int(C_ali4[3]), int(index_atom), int(index_atom+3)])
                #
                #if j+1==num_int_monomers_PS: ##penultimate monomer has different dihedral
                #   backbone_dihedral.append(int(C_ali4[3]))
                #   ##dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+15)])
                #   #dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+4)])
                #   #dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+5)])
                #else:  
                #   dihedral.append([1, int(C_ali4[3]), int(index_atom), int(index_atom+3), int(index_atom+16)])
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
        if(k == nB-1):
            mpb = int(lenp(t1[i*nB+k-1],t1[i*nB+k])/bondm)
        else:
            mpb = int(lenp(t1[i*nB+k+1],t1[i*nB+k])/bondm)
        #print "k", k
        #print "i",i
        for j in range(mpb):            
            if(j == mpb-1 and k == nB-1):    
                x = x + dx
                y = y + dy
                z = z + dz
                posc_ali4 = np.asarray([x,y,z])
                C_ali4 = np.asarray([posc_ali4[0], posc_ali4[1], posc_ali4[2], index_atom])
                atom.append([i, 13, 0.0031, C_ali4])
                index_atom +=1
                
                posh_ali1 = getposh_ali2(posc_ali4, vecv)
                H_ali1 = np.asarray([posh_ali1[0], posh_ali1[1],posh_ali1[2], index_atom])
                atom.append([i, 14, 0.0, H_ali1])
                index_atom +=1
                
                posh_ali2 = getposh_ali3(posc_ali4,vecv)
                H_ali2 = np.asarray([posh_ali2[0], posh_ali2[1],posh_ali2[2], index_atom])
                atom.append([i, 14, 0.0, H_ali2])
                index_atom +=1
                
                posc_ali5 = getposc_ali2(posc_ali4,vecm)
                C_ali5 = np.asarray([posc_ali5[0], posc_ali5[1], posc_ali5[2], index_atom])
                atom.append([i, 15, 0.0189, C_ali5])
                index_atom +=1
                
                posc_ali7 = getposc_aro5(posc_ali5, -vecmr)
                C_ali7 = np.asarray([posc_ali7[0], posc_ali7[1], posc_ali7[2], index_atom])
                atom.append([i, 13, 0.0031, C_ali7])
                index_atom +=1
                
                posh_ali3 = getposh_ali1(posc_ali7, -vec)
                H_ali3 = np.asarray([posh_ali3[0], posh_ali3[1], posh_ali3[2], index_atom])
                atom.append([i, 14, 0.0, H_ali3])
                index_atom +=1
                
                posh_ali7 = getposh_ali2(posc_ali7, vecv)
                H_ali7 = np.asarray([posh_ali7[0], posh_ali7[1], posh_ali7[2], index_atom])
                atom.append([i, 14, 0.0, H_ali7])
                index_atom +=1
                
                posh_ali8 = getposh_ali3(posc_ali7, vecv)
                H_ali8 = np.asarray([posh_ali8[0], posh_ali8[1], posh_ali8[2], index_atom])
                atom.append([i, 14, 0.0, H_ali8])
                index_atom +=1
                
                posc_ali6 = np.asarray([posc_ali5[0]+bondlength*vec[0],posc_ali5[1]+bondlength*vec[1],posc_ali5[2]+bondlength*vec[2]])
                C_ali6 = np.asarray([posc_ali6[0], posc_ali6[1], posc_ali6[2], index_atom])
                atom.append([i, 13, 0.0031, C_ali6])
                index_atom +=1
                
                vtem = np.asarray([posc_ali6[0]-posc_ali5[0],posc_ali6[1]-posc_ali5[1],posc_ali6[2]-posc_ali5[2]])
                vtem = vtem/lent(vtem)
                posh_ali4 = getposh_ali2(posc_ali6, vecm)
                H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                atom.append([i, 14, 0.0, H_ali4])
                index_atom +=1
                
                posh_ali5 = getposh_ali1(posc_ali6,vec)
                H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                atom.append([i, 14, 0.0, H_ali5])
                index_atom +=1
                
                posh_ali6 = getposh_ali1(posc_ali6,-vec)
                H_ali6 = np.asarray([posh_ali6[0], posh_ali6[1], posh_ali6[2], index_atom])
                atom.append([i, 14, 0.0, H_ali6])
                index_atom +=1
                
                posc_sid1 = getposc_aro6(posc_ali5,vecd)
                C_sid1 = np.asarray([posc_sid1[0], posc_sid1[1], posc_sid1[2], index_atom])
                atom.append([i, 23, 0.7474, C_sid1])
                index_atom +=1
                
                poso_sid1 = getposc_aro1(posc_sid1,vecm)
                O_sid1 = np.asarray([poso_sid1[0], poso_sid1[1], poso_sid1[2], index_atom])
                atom.append([i, 17, -0.5939, O_sid1])
                index_atom +=1
                
                poso_sid2 = getposc_aro5(posc_sid1,vecmr)
                O_sid2 = np.asarray([poso_sid2[0], poso_sid2[1], poso_sid2[2], index_atom])
                atom.append([i, 24, -0.4617, O_sid2])
                index_atom +=1
                
                posc_sid2 = getposc_aro4(poso_sid2,vecd)
                C_sid2 = np.asarray([posc_sid2[0],posc_sid2[1], posc_sid2[2], index_atom])
                atom.append([i, 25, 0.2800, C_sid2])
                index_atom +=1  
                
                posh_sid1 = getposh_ali1(posc_sid2,vec)
                H_sid1 = np.asarray([posh_sid1[0], posh_sid1[1], posh_sid1[2], index_atom])
                atom.append([i, 20, 0.0, H_sid1])
                index_atom +=1 
                
                posh_sid2 = getposh_ali2(posc_sid2,vecv)
                H_sid2 = np.asarray([posh_sid2[0], posh_sid2[1], posh_sid2[2], index_atom])
                atom.append([i, 20, 0.0, H_sid2])
                index_atom +=1
                
                posh_sid3 = getposh_ali3(posc_sid2,vecv)
                H_sid3 = np.asarray([posh_sid3[0], posh_sid3[1], posh_sid3[2], index_atom])
                atom.append([i, 20, 0.0, H_sid3])
                index_atom +=1 
                
                ## Bonding information
                bonding.append([13, int(C_ali4[3]), int(H_ali1[3])])
                bonding.append([13, int(C_ali4[3]), int(H_ali2[3])])            
                bonding.append([14, int(C_ali4[3]), int(C_ali5[3])])
                bonding.append([14, int(C_ali5[3]), int(C_ali6[3])])
                #bonding.append([1, int(C_ali5[3]), int(H_ali3[3])])
                bonding.append([14, int(C_ali5[3]), int(C_ali7[3])])
                bonding.append([13, int(C_ali7[3]), int(H_ali3[3])])
                bonding.append([13, int(C_ali7[3]), int(H_ali7[3])])
                bonding.append([13, int(C_ali7[3]), int(H_ali8[3])])
                bonding.append([13, int(C_ali6[3]), int(H_ali4[3])])
                bonding.append([13, int(C_ali6[3]), int(H_ali5[3])])
                bonding.append([13, int(C_ali6[3]), int(H_ali6[3])])
                bonding.append([15, int(C_ali5[3]), int(C_sid1[3])])
                bonding.append([16, int(C_sid1[3]), int(O_sid1[3])])
                bonding.append([17, int(C_sid1[3]), int(O_sid2[3])])
                bonding.append([18, int(O_sid2[3]), int(C_sid2[3])])
                bonding.append([19, int(C_sid2[3]), int(H_sid1[3])])
                bonding.append([19, int(C_sid2[3]), int(H_sid2[3])])
                bonding.append([19, int(C_sid2[3]), int(H_sid3[3])])
            
            ##Storing angle information of monomer
                angle.append([13, int(H_ali1[3]), int(C_ali4[3]), int(H_ali2[3])])
                angle.append([13, int(H_ali4[3]), int(C_ali6[3]), int(H_ali5[3])])
                angle.append([13, int(H_ali4[3]), int(C_ali6[3]), int(H_ali6[3])])
                angle.append([13, int(H_ali6[3]), int(C_ali6[3]), int(H_ali5[3])])
                angle.append([13, int(H_ali3[3]), int(C_ali7[3]), int(H_ali7[3])])
                angle.append([13, int(H_ali3[3]), int(C_ali7[3]), int(H_ali8[3])])
                angle.append([13, int(H_ali7[3]), int(C_ali7[3]), int(H_ali8[3])])
                angle.append([14, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3])])
                #angle.append([2, int(H_ali3[3]), int(C_ali5[3]), int(C_ali4[3])])
                #angle.append([2, int(H_ali3[3]), int(C_ali5[3]), int(C_ali6[3])])
                #angle.append([2, int(H_ali3[3]), int(C_ali5[3]), int(C_sid1[3])])
                
                angle.append([14, int(H_ali3[3]), int(C_ali7[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali7[3]), int(C_ali7[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali8[3]), int(C_ali7[3]), int(C_ali5[3])])
            
                angle.append([15, int(C_ali4[3]), int(C_ali5[3]), int(C_ali6[3])])
                angle.append([15, int(C_ali4[3]), int(C_ali5[3]), int(C_ali7[3])])
                angle.append([15, int(C_ali6[3]), int(C_ali5[3]), int(C_ali7[3])])
            
                angle.append([16, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3])])
                angle.append([16, int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3])])
            
                angle.append([16, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                angle.append([17, int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                angle.append([18, int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                angle.append([19, int(O_sid1[3]), int(C_sid1[3]), int(O_sid2[3])])
                angle.append([20, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                angle.append([21, int(O_sid2[3]), int(C_sid2[3]), int(H_sid1[3])])
                angle.append([21, int(O_sid2[3]), int(C_sid2[3]), int(H_sid2[3])])
                angle.append([21, int(O_sid2[3]), int(C_sid2[3]), int(H_sid3[3])])
                angle.append([22, int(H_sid1[3]), int(C_sid2[3]), int(H_sid2[3])])
                angle.append([22, int(H_sid1[3]), int(C_sid2[3]), int(H_sid3[3])])
                angle.append([22, int(H_sid2[3]), int(C_sid2[3]), int(H_sid3[3])])
                
            ##Storing dihedral information for last monomer     
                dihedral.append([3, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_ali7[3])])
                dihedral.append([3, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_ali7[3])])
                dihedral.append([3, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali7[3])])
                dihedral.append([3, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali7[3])])
                dihedral.append([3, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali7[3])])
                dihedral.append([3, int(H_ali3[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali7[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali8[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali3[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(H_ali7[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(H_ali8[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_ali6[3])])    
            
                
                dihedral.append([4, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([4, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([4, int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([5, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([5, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([5, int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
            
                
                dihedral.append([6, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([6, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([6, int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                
                dihedral.append([7, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([7, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([7, int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                
                dihedral.append([9, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali3[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali7[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali8[3]), int(C_ali7[3]), int(C_ali5[3]), int(C_sid1[3])])
                
                dihedral.append([10, int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                dihedral.append([11, int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
            
                
                dihedral.append([13, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3]), int(H_sid1[3])])
                dihedral.append([13, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3]), int(H_sid2[3])])
                dihedral.append([13, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3]), int(H_sid3[3])])
                
                dihedral.append([14, int(O_sid1[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                dihedral.append([15, int(O_sid1[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
            
                dihedral.append([16, int(C_sid1[3]), int(O_sid1[3]), int(C_sid2[3]), int(O_sid2[3])])
            
                #Adding backbone dihedrals
                backbone_dihedral = np.asarray(backbone_dihedral)
                for i in range(len(backbone_dihedral)-4):
                    dihedral.append([1, backbone_dihedral[i], backbone_dihedral[i+1], backbone_dihedral[i+2], backbone_dihedral[i+3]])        
            else:
                if(j==0):
                    #print "for internal monomer,j=0", i*nB+k+1
                   # print " for internal monomer, j = 0", j ,k,i
                    if(k == nB-1):
                        vec = np.asarray([t1[i*nB+k][0]-t1[i*nB+k-1][0],t1[i*nB+k][1]-t1[i*nB+k-1][1],t1[i*nB+k][2]-t1[i*nB+k-1][2]])
                    else:
                        vec = np.asarray([t1[i*nB+k+1][0]-t1[i*nB+k][0],t1[i*nB+k+1][1]-t1[i*nB+k][1],t1[i*nB+k+1][2]-t1[i*nB+k][2]])
                    lentv = lent(vec)
                    vec = vec/lent(vec)
                    #a unit vector which is vertical t vec
                    vecv = dirv(vec)
                    vecd = dirvd(vec, vecv)
                    vecm = 1.732*vec + vecd
                    vecm = vecm/lent(vecm)
                    vecmr = vecd - 1.732*vec
                    vecmr = vecmr/lent(vecmr)
                    dx = vec[0]*lentv/mpb
                    dy = vec[1]*lentv/mpb
                    dz = vec[2]*lentv/mpb           
                    x = t1[i*nB+k][0]
                    y = t1[i*nB+k][1]
                    z = t1[i*nB+k][2]
                else:
                   # print " for internal monomer", j ,k, i*nB+k
                    x = x + dx
                    y = y + dy
                    z = z + dz
                    
                posc_ali4 = np.asarray([x,y,z])
                C_ali4 = np.asarray([posc_ali4[0], posc_ali4[1], posc_ali4[2], index_atom])
                atom.append([i, 21, 0.0051, C_ali4])
                index_atom +=1
                
                posh_ali1 = getposh_ali2(posc_ali4, vecv)
                H_ali1 = np.asarray([posh_ali1[0], posh_ali1[1],posh_ali1[2], index_atom])
                atom.append([i, 14, 0.0, H_ali1])
                index_atom +=1
                
                posh_ali2 = getposh_ali3(posc_ali4,vecv)
                H_ali2 = np.asarray([posh_ali2[0], posh_ali2[1],posh_ali2[2], index_atom])
                atom.append([i, 14, 0.0, H_ali2])
                index_atom +=1
                
                posc_ali5 = getposc_ali2(posc_ali4,vecm)
                C_ali5 = np.asarray([posc_ali5[0], posc_ali5[1], posc_ali5[2], index_atom])
                atom.append([i, 15, 0.0189, C_ali5])
                index_atom +=1
                C_ali2 = C_ali5
                
                posc_ali6 = getposh_ali4(posc_ali5,vecv)
                C_ali6 = np.asarray([posc_ali6[0], posc_ali6[1], posc_ali6[2], index_atom])
                atom.append([i, 21, 0.0051, C_ali6])
                index_atom +=1
                
                vtem = np.asarray([posc_ali6[0]-posc_ali5[0],posc_ali6[1]-posc_ali5[1],posc_ali6[2]-posc_ali5[2]])
                vtem = vtem/lent(vtem)
                posh_ali4 = getposh_ali2(posc_ali6, vtem)
                H_ali4 = np.asarray([posh_ali4[0], posh_ali4[1], posh_ali4[2], index_atom])
                atom.append([i, 14, 0.0, H_ali4])
                index_atom +=1
                
                posh_ali5 = getposh_ali1(posc_ali6,vec)
                H_ali5 = np.asarray([posh_ali5[0], posh_ali5[1], posh_ali5[2], index_atom])
                atom.append([i, 14, 0.0, H_ali5])
                index_atom +=1
                
                posh_ali6 = getposh_ali1(posc_ali6,-vec)
                H_ali6 = np.asarray([posh_ali6[0], posh_ali6[1], posh_ali6[2], index_atom])
                atom.append([i, 14, 0.0, H_ali6])
                index_atom +=1
                
                posc_sid1 = getposc_aro6(posc_ali5,vecd)
                C_sid1 = np.asarray([posc_sid1[0], posc_sid1[1], posc_sid1[2], index_atom])
                atom.append([i, 16, 0.7464, C_sid1])
                index_atom +=1
                
                poso_sid1 = getposc_aro1(posc_sid1,vecm)
                O_sid1 = np.asarray([poso_sid1[0], poso_sid1[1], poso_sid1[2], index_atom])
                atom.append([i, 17, -0.5939, O_sid1])
                index_atom +=1
                
                poso_sid2 = getposc_aro5(posc_sid1,vecmr)
                O_sid2 = np.asarray([poso_sid2[0], poso_sid2[1], poso_sid2[2], index_atom])
                atom.append([i, 18, -0.4617, O_sid2])
                index_atom +=1
                
                posc_sid2 = getposc_aro4(poso_sid2,vecd)
                C_sid2 = np.asarray([posc_sid2[0],posc_sid2[1], posc_sid2[2], index_atom])
                atom.append([i, 22, 0.2801, C_sid2])
                index_atom +=1  
                
                posh_sid1 = getposh_ali1(posc_sid2,vec)
                H_sid1 = np.asarray([posh_sid1[0], posh_sid1[1], posh_sid1[2], index_atom])
                atom.append([i, 20, 0.0, H_sid1])
                index_atom +=1 
                
                posh_sid2 = getposh_ali2(posc_sid2,vecv)
                H_sid2 = np.asarray([posh_sid2[0], posh_sid2[1], posh_sid2[2], index_atom])
                atom.append([i, 20, 0.0, H_sid2])
                index_atom +=1
                
                posh_sid3 = getposh_ali3(posc_sid2,vecv)
                H_sid3 = np.asarray([posh_sid3[0], posh_sid3[1], posh_sid3[2], index_atom])
                atom.append([i, 20, 0.0, H_sid3])
                index_atom +=1 
                
                ### Bonding information for monomer
                bonding.append([13, int(C_ali4[3]), int(H_ali1[3])])
                bonding.append([13, int(C_ali4[3]), int(H_ali2[3])])            
                bonding.append([14, int(C_ali4[3]), int(C_ali5[3])])
                bonding.append([14, int(C_ali5[3]), int(C_ali6[3])])
                bonding.append([14, int(C_ali5[3]), int(index_atom)]) ##bonding to next monomer
                bonding.append([13, int(C_ali6[3]), int(H_ali4[3])])
                bonding.append([13, int(C_ali6[3]), int(H_ali5[3])])
                bonding.append([13, int(C_ali6[3]), int(H_ali6[3])])
                bonding.append([15, int(C_ali5[3]), int(C_sid1[3])])
                bonding.append([16, int(C_sid1[3]), int(O_sid1[3])])
                bonding.append([17, int(C_sid1[3]), int(O_sid2[3])])
                bonding.append([18, int(O_sid2[3]), int(C_sid2[3])])
                bonding.append([19, int(C_sid2[3]), int(H_sid1[3])])
                bonding.append([19, int(C_sid2[3]), int(H_sid2[3])])
                bonding.append([19, int(C_sid2[3]), int(H_sid3[3])])
                
                
        ##Storing angle information of monomer
                angle.append([13, int(H_ali1[3]), int(C_ali4[3]), int(H_ali2[3])])
                angle.append([13, int(H_ali4[3]), int(C_ali6[3]), int(H_ali5[3])])
                angle.append([13, int(H_ali4[3]), int(C_ali6[3]), int(H_ali6[3])])
                angle.append([13, int(H_ali6[3]), int(C_ali6[3]), int(H_ali5[3])])
                angle.append([14, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3])])
                angle.append([14, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3])])
                angle.append([14, int(index_atom+1), int(index_atom), int(C_ali5[3])])
                angle.append([14, int(index_atom+2), int(index_atom), int(C_ali5[3])])
                angle.append([15, int(C_ali4[3]), int(C_ali5[3]), int(C_ali6[3])])
                angle.append([15, int(C_ali4[3]), int(C_ali5[3]), int(index_atom)])
                angle.append([15, int(C_ali5[3]), int(index_atom), int(index_atom+3)])
                angle.append([15, int(C_ali6[3]), int(C_ali5[3]), int(index_atom)])
                angle.append([16, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3])])
                angle.append([16, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                angle.append([16, int(index_atom), int(C_ali5[3]), int(C_sid1[3])])
                angle.append([17, int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                angle.append([18, int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                angle.append([19, int(O_sid1[3]), int(C_sid1[3]), int(O_sid2[3])])
                angle.append([20, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                angle.append([21, int(O_sid2[3]), int(C_sid2[3]), int(H_sid1[3])])
                angle.append([21, int(O_sid2[3]), int(C_sid2[3]), int(H_sid2[3])])
                angle.append([21, int(O_sid2[3]), int(C_sid2[3]), int(H_sid3[3])])
                angle.append([22, int(H_sid1[3]), int(C_sid2[3]), int(H_sid2[3])])
                angle.append([22, int(H_sid1[3]), int(C_sid2[3]), int(H_sid3[3])])
                angle.append([22, int(H_sid2[3]), int(C_sid2[3]), int(H_sid3[3])])
                
            ##Storing dihedral information of monomer
                dihedral.append([3, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3]), int(index_atom)])
                dihedral.append([3, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3]), int(index_atom)])
                dihedral.append([3, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3]), int(index_atom)])
                dihedral.append([3, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3]), int(index_atom)])
                dihedral.append([3, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3]), int(index_atom)])
                dihedral.append([3, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(index_atom+1), int(index_atom), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(index_atom+1), int(index_atom), int(C_ali5[3]), int(C_ali6[3])])
                dihedral.append([3, int(index_atom+2), int(index_atom), int(C_ali5[3]), int(C_ali4[3])])
                dihedral.append([3, int(index_atom+2), int(index_atom), int(C_ali5[3]), int(C_ali6[3])])
                
                dihedral.append([4, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([4, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([4, int(index_atom), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([5, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([5, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                dihedral.append([5, int(index_atom), int(C_ali5[3]), int(C_sid1[3]), int(O_sid1[3])])
                
                dihedral.append([6, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([6, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([6, int(index_atom), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                
                dihedral.append([7, int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([7, int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                dihedral.append([7, int(index_atom), int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3])])
                
                dihedral.append([8, int(C_ali4[3]), int(C_ali5[3]), int(index_atom), int(index_atom+3)])
                dihedral.append([8, int(C_ali6[3]), int(C_ali5[3]), int(index_atom), int(index_atom+3)])
                dihedral.append([8, int(C_ali5[3]), int(index_atom), int(index_atom+3), int(index_atom+4)]) 
                if j + 1< mpb or k +1 < nB :
                    dihedral.append([8, int(C_ali5[3]), int(index_atom), int(index_atom+3), int(index_atom+15)]) ###Check for last index after monomer is set up
                    dihedral.append([12, int(C_ali5[3]), int(index_atom), int(index_atom+3), int(index_atom+8)])
                else:
                    #dihedral.append([6, int(C_ali5[3]), int(index_atom), int(index_atom+3), int(index_atom+4)]) ###Check for last index after monomer is set up
                    dihedral.append([8, int(C_ali5[3]), int(index_atom), int(index_atom+3), int(index_atom+8)]) ###Check for last index after monomer is set up
                    dihedral.append([12, int(C_ali5[3]), int(index_atom), int(index_atom+3), int(index_atom+12)])
                
                dihedral.append([9, int(H_ali1[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali2[3]), int(C_ali4[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali4[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali5[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(H_ali6[3]), int(C_ali6[3]), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(index_atom+1), int(index_atom), int(C_ali5[3]), int(C_sid1[3])])
                dihedral.append([9, int(index_atom+2), int(index_atom), int(C_ali5[3]), int(C_sid1[3])])
                
                dihedral.append([10, int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                dihedral.append([11, int(C_ali5[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                
                dihedral.append([12, int(index_atom+3), int(index_atom), int(C_ali5[3]), int(C_sid1[3])])
                
                
                dihedral.append([13, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3]), int(H_sid1[3])])
                dihedral.append([13, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3]), int(H_sid2[3])])
                dihedral.append([13, int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3]), int(H_sid3[3])])
                
                dihedral.append([14, int(O_sid1[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
                dihedral.append([15, int(O_sid1[3]), int(C_sid1[3]), int(O_sid2[3]), int(C_sid2[3])])
        
                dihedral.append([16, int(C_sid1[3]), int(O_sid1[3]), int(C_sid2[3]), int(O_sid2[3])])


total_atoms = len(atom)
total_bonds = len(bonding)
total_angles = len(angle)
total_dihedrals = len(dihedral)
total_impropers = len(improper)
#total_atoms = num_chains*(17 + 16*num_int_monomers_PS + 15*num_int_monomers_P2VP + 16)
#total_bonds = num_chains*(18 + 17*num_int_monomers_PS + 16*num_int_monomers_P2VP + 16)
#total_angles = num_chains*(33 + 30*num_int_monomers_PS + 28*num_int_monomers_P2VP + 25)
#total_dihedrals = num_chains*(5 + 2*(num_int_monomers_PS - 1) + 2 + 2*(num_int_monomers_P2VP - 1) + 3 + 0)
#total_impropers = num_chains*(12 + 12*num_int_monomers_PS + 0)

f = open("PS_PMMA_50.data","w")
###Headers
with open("PS_PMMA_50.data", "a") as out_file:
    out_file.write('PS with ' + str(num_chains) +' chains and ' + str(num_int_monomers_PS) + ' internal PS monomers per chain and ' + str(num_int_monomers_PMMA) +' int PMMA monomers per chain'+ 'and ' + str(number_ions) + ' ions' +'\n'
    + '\n' +
    str(total_atoms) + ' ' + 'atoms' + '\n'
    + str(total_bonds) + ' ' + 'bonds' + '\n'
    + str(total_angles) + ' ' + 'angles' + '\n'
    + str(total_dihedrals) + ' ' + 'dihedrals' + '\n'
    + str(total_impropers) + ' ' + 'impropers' + '\n'
    + '\n'
    + """26 atom types
19 bond types
22 angle types
16 dihedral types
2 improper types

0 250 xlo xhi
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
11 32.065
12 15.99 
13 12.0
14 1.00782
15 12.0
16 12.0
17 15.9949
18 15.9949
19 12.0
20 1.00782
21 12.0
22 12.0
23 12.0
24 15.9949
25 12.0
26 22.99

Bond Coeffs

1 310.44 1.53 #134 1.53
2 317.45 1.51 #158.5 1.51
3 469.47 1.39 #234.5 1.39
4 239.00 1.1  #170 1.1
5 239.00 1.08 #183.5 1.08
6  310.44 1.53 #134 1.53
7  317.45 1.51 #158.5 1.51
8  469.47 1.39 #234.5 1.39
9  239.00 1.1  #170 1.1
10 239.00 1.08 #183.5 1.08
11 167.305 1.85
12 270.08 1.66 
13 360.24091778202677 1.095 
14 310.2074569789675 1.52 
15 330.2208413001912 1.5 
16 999.6684990439771 1.204 
17 411.5752868068834 1.343 
18 320.2141491395794 1.41 
19 410.274378585086 1.092
 
Angle Coeffs

1 43.8456 109.45
2 36.6157 109.45
3 57.6362 109.45
4 45.0048 120
5 50.0478 120
6  43.8456 109.45
7  36.6157 109.45
8  57.6362 109.45
9  45.0048 120
10 50.0478 120
11 28.08   103.8
12 50.31   114.8
13 40.02676864244742 109.5 
14 45.030114722753346 112.6 
15 45.030114722753346 113.5 
16 40.02676864244742 111.5 
17 59.039483747609935 125.4 
18 50.03346080305927 111 
19 103.0689292543021 122.5 
20 50.633867112810705 114 
21 60.040152963671126 110 
22 45.030114722753346 109.5 
23 74.34 120

Dihedral Coeffs

1 harmonic 1.43403 -1 3
2 opls 0.741 -2.25 1.39 0
3  charmm 0.050033460803059274   3 180 0 
4  charmm -0.20513742829827916   2 90  0  
5  charmm -0.14509679732313574   3 180 0 
6  charmm -0.22515057361376672   2 90  0  
7  charmm -0.16510994263862333   3 180 0 
8  charmm 0.22515057361376672    1 180 0 
9  charmm -0.050033460803059274  3 180 0 
10 charmm -0.625418260038241     1 180 0 
11 charmm -1.2758484703632886    2  90 0  
12 charmm  0.0750501912045889    2  90 0 
13 charmm  -0.0750501912045889   3 180 0
14 charmm -2.4066085086042066    2  90 0 
15 charmm  -0.5003346080305927   3 180 0
16 charmm  2.001338432122371     2   0 0 

Improper Coeffs

1 20.0048 0
2 47.8 0""" + '\n' + '\n' + 'Atoms' +'\n' +'\n')

## add Atoms

for i in range(len(atom)):
    write_data_atom(i+1, atom[i])
#with open("PS_P2VP_sac.data", "r") as in_file:
#   with open("PS_P2VP_50.data", "a") as out_file:
#       out_file.write('\n' + "Atoms" + '\n' + '\n')
#       all_coords = in_file.readlines()
#       out_file.write(all_coords)


## Adding Bonding information
bonding=np.asarray(bonding) 
with open("PS_PMMA_50.data", "a") as out_file:
    out_file.write('\n'+"Bonds" + '\n' + '\n')
for i in range(len(bonding)):
    write_data_bond(i+1, bonding[i])
## Adding Angle information
angle = np.asarray(angle)
with open("PS_PMMA_50.data", "a") as out_file:
    out_file.write('\n' + "Angles" + '\n' + '\n')
for i in range(len(angle)):
    write_data_angle(i+1, angle[i])
##Adding dihedral information
dihedral = np.asarray(dihedral)
with open("PS_PMMA_50.data","a") as out_file:
    out_file.write('\n' + "Dihedrals" + '\n' + '\n')
for i in range(len(dihedral)):
    write_data_dihedral(i+1, dihedral[i])
##Adding improper information
improper = np.asarray(improper)
with open("PS_PMMA_50.data", "a") as out_file:
    out_file.write('\n' + "Impropers" + '\n' + '\n')
for i in range(len(improper)):
    write_data_improper(i+1, improper[i])
                                        
                                
                    
