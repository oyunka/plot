# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys
import math
from matplotlib import pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import lagrange
argv = sys.argv
argc = len(argv)

if(argc != 2):
    print("one input is required.")
    exit()

fname = argv[1]#"H_ecs_r40toinf_th15_s2t20.surff"

pi = np.pi
data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None, skiprows=3)

fname = "betha_parameters_each_k.txt"
betha_data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None)
wang_data = pd.read_csv("wang.txt", comment="#", delim_whitespace=True, header=None, skiprows=1)


num_krad = 400
num_kang = 24

krad = data.ix[:,0]
kang = data.ix[:,1]
kval = data.ix[:,2]

betha1 = betha_data.ix[:,0]
betha2 = betha_data.ix[:,1]
betha3 = betha_data.ix[:,2]
betha4 = betha_data.ix[:,3]

cos_ang = wang_data.ix[:,0]
kang_rad = wang_data.ix[:,1]
weight = wang_data.ix[:,2]

erad= krad * krad * 27.2 * 0.5
prob = krad* krad* abs(kval) / (krad)



krad = np.reshape(krad, (num_krad, num_kang))
kang = np.reshape(kang, (num_krad, num_kang))



cmap = cm.jet

num_krad_half = 200
orbital_2s = 0


def Legendre1():
    return  cos_ang

def Legendre2():
    first = np.zeros((num_kang), dtype = 'float')
    for m in range(0, num_kang):
        first[m] += 0.5 * (3*(cos_ang[m])**2 - 1)   
    return  first

def Legendre3():
    first = np.zeros((num_kang), dtype = 'float')

    for m in range(0, num_kang):
        first[m] +=0.5 * (5 * pow(cos_ang[m], 3)  - 3 * cos_ang[m])
    return  first

def Legendre4():
    first = np.zeros((num_kang), dtype = 'float')

    for m in range(0, num_kang):
        first[m] += 1/8 * (35 * pow(cos_ang[m], 4) - 30 * pow (cos_ang[m], 2) + 3)
    return  first


I = np.empty((num_krad, num_kang), dtype='float')



legendre1 = Legendre1()
legendre2 = Legendre2()
legendre3 = Legendre3()
legendre4 = Legendre4()

for i in range(0, num_krad):
    for j in range(0, num_kang):
        I[i, j] = betha1[i] * legendre1[j]+ betha2[i] * legendre2[j] + betha3[i] * legendre3[j] + betha4[i] * legendre4[j]

###############################################################
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




# fig = plt.figure(figsize=(8, 6))

# ###############################################################
# if(1):
#     ax1 = fig.add_subplot(111, projection='polar')
#     ax1.set_theta_zero_location('N')
#     ax1.set_theta_direction(-1)     
    
#     cmin = 0
#     cmax = 0.1
    
   
#     im3 = ax1.pcolormesh(kang, erad, prob, cmap=cmap,  shading='gouraud', rasterized=True)
#     im3.set_clim(cmin,cmax)
#     im4 = ax1.pcolormesh(-kang, erad, prob, cmap=cmap,  shading='gouraud', rasterized=True)
#     im4.set_clim(cmin, cmax)
    
#     ax1.axis("image")
#     plt.colorbar(im3, ax=ax1)
#     plt.grid(lw=0.4, color="black", linestyle="-")
    
#     ax1.set_rmin(0)
#     ax1.set_rmax(40)
#     ax1.set_title("Angular distribution of energy spectrum (linear scale) \n\n")
# plt.show()

fig = plt.figure(figsize=(16, 6))

###############################################################
if(1):
    ax1 = fig.add_subplot(121, projection='polar')
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)     
    
    cmin = 0
    cmax = 0.03
    
    im3 = ax1.pcolormesh(kang, erad, prob, cmap=cmap,  shading='gouraud', rasterized=True)
    im3.set_clim(cmin,cmax)
    im4 = ax1.pcolormesh(-kang, erad, prob, cmap=cmap,  shading='gouraud', rasterized=True)
    im4.set_clim(cmin, cmax)
    
    ax1.axis("image")
    plt.colorbar(im3, ax=ax1)
    plt.grid(lw=0.4, color="black", linestyle="-")
    
    ax1.set_rmin(0)
    ax1.set_rmax(8)
    ax1.set_title("angular distribution of energy spectrum (linear scale) \n\n")


###############################################################
ax2 = fig.add_subplot(122, projection='polar')
ax2.set_theta_zero_location('N')
ax2.set_theta_direction(-1)

cmin = 10**(-8.0)
cmax = 10**(10.0)

im = ax2.pcolormesh(kang, erad, prob, cmap=cmap, norm=LogNorm(), shading='gouraud', rasterized=True)
im.set_clim(cmin,cmax)
im2 = ax2.pcolormesh(-kang, erad, prob, cmap=cmap, norm=LogNorm(), shading='gouraud', rasterized=True)
im2.set_clim(cmin, cmax)

ax2.axis("image")
cax2 = fig.add_axes([0.92, 0.1, 0.015, 0.8])   
cb2 = plt.colorbar(im, cax2, ax=ax2)
#plt.colorbar(im, ax=ax2)

cb2.set_label("yield (arb. units)", size=16)
cb2.ax.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)

#ax2.set_xticklabels(["a", "b"])
ax2.grid(lw=0.4, color="black", linestyle="-")
#ax2.set_yticklabels(["50 eV", "100 eV", "150 eV", "200 eV"])



ax2.set_rmin(0)
ax2.set_rmax(40)
ax2.set_title("angular distribution of energy spectrum (log scale) \n\n ")

###############################################################


if(0):
    cmin = -15
    cmax = -4
    im = ax.pcolormesh(kang, erad, np.log(prob), cmap=cmap)
    im.set_clim(cmin,cmax)
    im2 = ax.pcolormesh(-kang, erad, np.log(prob), cmap=cmap)
    im2.set_clim(cmin, cmax)
plt.colorbar(im2)
plt.axis('equal')
plt.savefig(fname+'.pdf', bbox_inches='tight')
# plt.show()