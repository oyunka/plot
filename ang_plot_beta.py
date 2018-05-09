# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import sys
import math
from scipy.interpolate import lagrange
argv = sys.argv
argc = len(argv)

if(argc != 2):
    print("one input is required.")
    exit()

fname = argv[1]#"H_ecs_r40toinf_th15_s2t20.surff"

pi = np.pi

def generate_cmap(colors):
    """自分で定義したカラーマップを返す"""
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


with open(fname, "r") as fp:
    fp.readline()
    num_krad = int(fp.readline().strip("\n"))
    num_kang = int(fp.readline().strip("\n"))

data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None, skiprows=3)
wang_data = pd.read_csv("wang.txt", comment="#", delim_whitespace=True, header=None, skiprows=1)



krad = data.ix[:,0]
kang = data.ix[:,1]
kval = data.ix[:,2]

cos_ang = wang_data.ix[:,0]
kang_rad = wang_data.ix[:,1]
weight = wang_data.ix[:,2]


def_krad = krad
def_kang = kang

#kval = kr + 1.0j * ki

#R, Theta = np.meshgrid(r, theta) 

krad = np.reshape(krad, (num_krad, num_kang))
kang = np.reshape(kang, (num_krad, num_kang))
kval = np.reshape(kval, (num_krad, num_kang))


print(weight.size)

if(0):
    #new_kcosang = np.linspace(-1, 1, 20)
    new_kcosang = np.cos(kang[0,:])
    new_krad = np.empty((krad.shape[0], new_kcosang.size), dtype='float')
    new_kang = np.empty((krad.shape[0], new_kcosang.size), dtype='float')
    new_kval = np.empty((krad.shape[0], new_kcosang.size), dtype='complex')
    
    def angle_interpolate2D(krad, kang, kval, new_kcosang, new_krad, new_kang, new_kval):
        print("in", krad.shape[0])
        for ikrad in range(krad.shape[0]):
            print(ikrad, "done", end=" ")
            new_krad[ikrad, :] = krad[ikrad][0] * np.ones(new_kcosang.size)
            new_kang[ikrad, :] = np.arccos(new_kcosang[:])
            linter = lagrange(np.cos(kang[ikrad, :]), kval[ikrad, :])
            new_kval[ikrad, :] = linter(new_kcosang[:])
            if(1):
                print(np.cos(kang[ikrad, :]) - new_kcosang[:])
                print(kval[ikrad, :] - new_kval[ikrad, :])

        return new_krad, new_kang, new_kval

    (krad, kang, kval) = angle_interpolate2D(krad, kang, kval, new_kcosang, new_krad, new_kang, new_kval)
    print("conversion done")

# new_kcosang = np.cos(kang[0,:]) # ラジアンpiからpiまで入力


# print(kang)
# print(.
# kang.size)
# print(kval.size)
# print(krad)
# print(kang)
# print(kval)

## DEBUGING ------
#thresh =  2.8429890500E+00
#kang = kang * (kang > thresh) - 0 * (kang <= thresh)
#print(kang)
## ---------------
#krad = np.reshape(krad, (num_kang, num_krad))
#kang = np.reshape(kang, (num_kang, num_krad))
#kval = np.reshape(kval, (num_kang, num_krad))

#krad = (krad > 10**-8) * krad + (krad <= 10**-8) * 10**-8 
# erad = np.empty((krad.shape[0], num_kang), dtype='float')
# prob = np.empty((krad.shape[0], num_kang), dtype='float')

X1= (krad * krad * 27.2 * 0.5) * np.cos(kang)
X2 = (krad * krad * 27.2 * 0.5) * np.sin(kang)
erad= krad * krad * 27.2 * 0.5
prob = krad* krad* abs(kval) / (krad)



#######################Left Right Asymmetry############################

num_krad_half = 200
momentum = krad[:,0]
orbital_2s = 0
dk = np.zeros((num_krad), dtype = 'float')
mval = np.zeros((num_krad_half,  num_kang),dtype = 'float')

for i in range(0,num_krad):
    if(i < num_krad):
        dk[i] = momentum[i] - momentum[i-1]   
    if(i == 0) :
        dk[i] = momentum[i] - 0.0


for i in range(orbital_2s, num_krad_half):
    for k in range(0, num_kang):
        mval[i,k]  =   kval[i, k] * dk[i] * momentum[i] * momentum[i]
   
right = left  = 0

for l in range(orbital_2s,num_krad_half):
    for m in range(0, 12):
        left += mval[l, m]* weight[m]
    for m in range(12, 24):
        right += mval[l, m]* weight[m]

lr_asymmetry = (left - right)/(right + left)
print ('Asymmetry is',str(lr_asymmetry))
f = open('asymmetry.txt','a')
f.write(str(lr_asymmetry))
f.write('\n')
f.close()



#######################Legendre Polynomial fitting2 ############################
# print(cos_ang)

def Legendre0():
    first = np.zeros((num_krad_half), dtype = 'float')

    for l in range(orbital_2s, num_krad_half):
        for m in range(0, num_kang):
            first[l] += kval[l , m] * weight[m]
    all = 0

    for i in range(orbital_2s, num_krad_half):
        all +=   first[i]  * dk[i] * momentum[i] * momentum[i]
    return  all / 2


def Legendre1():
    first = np.zeros((num_krad_half), dtype = 'float')

    for l in range(orbital_2s,num_krad_half):
        for m in range(0, num_kang):
            first[l] += kval[l , m] * weight[m] * cos_ang[m]
    all = 0

    for i in range(orbital_2s, num_krad_half):
        all +=   first[i]  * dk[i] * momentum[i] * momentum[i]
    return  all * 3/2

def Legendre2():
    first = np.zeros((num_krad_half), dtype = 'float')

    for l in range(orbital_2s,num_krad_half):
        for m in range(0, num_kang):
            first[l] += kval[l , m] * weight[m] * 0.5 * (3*(cos_ang[m])**2 - 1)
    all = 0

    for i in range(orbital_2s, num_krad_half):
        all +=   first[i]  * dk[i] * momentum[i] * momentum[i]
# * dk[i] * momentum[i] * momentum[i]
    
    return  all * 5/2

def Legendre3():
    first = np.zeros((num_krad_half), dtype = 'float')

    for l in range(orbital_2s,num_krad_half):
        for m in range(0, num_kang):
            first[l] += kval[l , m] * weight[m] * 0.5 * (5* pow(cos_ang[m], 3) - 3 * cos_ang[m])
    all = 0
    for i in range(orbital_2s, num_krad_half):
        all +=   first[i]  * dk[i] * momentum[i] * momentum[i]
# * dk[i] * momentum[i] * momentum[i]


    return all * 7/2

def Legendre4():
    first = np.zeros((num_krad_half), dtype = 'float')

    for l in range(orbital_2s,num_krad_half):
        for m in range(0, num_kang):
            first[l] += kval[l , m] * weight[m] * 1/8 * (35 * pow(cos_ang[m], 4) - 30 * pow (cos_ang[m], 2) + 3)
    all = 0

    for i in range(orbital_2s, num_krad_half):
        all +=   first[i]   * dk[i] * momentum[i] * momentum[i]
# * dk[i] * momentum[i] * momentum[i]
  
    return all * 9/2
   

betha0 = Legendre0()
betha1 = Legendre1()/betha0
betha2 = Legendre2()/betha0
betha3 = Legendre3()/betha0
betha4 = Legendre4()/betha0

print(betha0,Legendre1(),betha1)

f = open('betha_parameters.txt','w')

f.write(str(betha1))
f.write('\t')
f.write(str(betha2))
f.write('\t')
f.write(str(betha3))
f.write('\t')
f.write(str(betha4))
f.write('\n')
f.close()


