# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import sys

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
w = np.fromfile("weight.txt", dtype=float, count=-1, sep=" ")
# w48 = pd.read_csv("weight48.txt", comment="#", delim_whitespace=True, header=None)

krad = data.ix[:,0]
kang = data.ix[:,1]
kval = data.ix[:,2]
#ki = data.ix[:,3]

def_krad = krad
def_kang = kang



# w = w48.ix[:,2]
weight = np.reshape(w, (num_kang, 1))
#R, Theta = np.meshgrid(r, theta)
krad = np.reshape(krad, (num_krad, num_kang))
kang = np.reshape(kang, (num_krad, num_kang))
kval = np.reshape(kval, (num_krad, num_kang))

# for i in range(weight.shape[0]):
#     for ikval in range(kval.shape[0]):
#          kval[ikval, i] = kval[ikval, i]* weight[i]

#dummy = np.ones((krad.shape[0],1), dtype="float") * np.pi * 1.01
#print(dummy.shape)
#krad = np.c_[dummy, krad]
#kang = np.c_[dummy, kang]
#kval = np.c_[dummy, kval]

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


# print(X1.size)
#cmap = cm.spectral
#cmap = cm.gist_rainbow
#cmap = cm.gist_ncar
#cmap = cm.nipy_spectral
#cmap = cm.gnuplot
#cmap = cm.gnuplot2
#cmap = cm.Paired # <- best
#cmap = cm.Spectral
#cmap = cm.Set1
#cmap = generate_cmap(['blue', 'lightblue', "yellow", "orange", 'red'])
cmap = cm.jet

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
plt.colorbar(im, ax=ax2)

cb2.set_label("yield (arb. units)", size=16)
cb2.ax.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)

ax2.set_xticklabels(["a", "b"])
ax2.grid(lw=0.4, color="black", linestyle="-")
ax2.set_yticklabels(["50 eV", "100 eV", "150 eV", "200 eV"])



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

num_krad_half = int(num_krad/2)
num_kang_half = int(num_kang/2)
print(num_krad)
seperation = 0
print(num_krad_half, num_kang_half)
momentum = krad[:,0]

dk = np.zeros((num_krad), dtype = 'float')
weighted_kval = np.zeros((num_krad, num_kang), dtype = 'float')



pie1 = pie2 = pie3 = pie4 = 0
for i in range(0, num_krad):
    if(i < num_krad):
        dk[i] = momentum[i] - momentum[i-1]
    if(i == 0) :
        dk[i] = momentum[i] - 0.0

for i in range(0, num_krad):
    for j in range(0, num_kang):
        weighted_kval[i, j] = kval[i, j] * weight[j]

print(kval.shape)


for i in range(0, num_krad_half):
    for j in range(0, num_kang_half):
        pie1 +=  weighted_kval[i, j] * momentum[i] * momentum[i] * dk[i]
    for j in range(num_kang_half, num_kang):
        pie2 +=  weighted_kval[i, j] * momentum[i] * momentum[i] * dk[i]



for i in range(num_kang_half, num_krad):
    for j in range(0, num_kang_half):
        pie3 +=  weighted_kval[i, j] * momentum[i] * momentum[i] * dk[i]
    for j in range(num_kang_half, num_kang):
        pie4 +=  weighted_kval[i, j] * momentum[i] * momentum[i] * dk[i]

print(pie1, pie2)



lr_asymmetry1 = (pie1 - pie2)/(pie1 + pie2)
lr_asymmetry2 = (pie3 - pie4)/(pie3 + pie4)
print(lr_asymmetry1, lr_asymmetry2)


print(weighted_kval[57,: ],kang[57, :] )
y = weighted_kval[57, :]
x = kang[57, :]
plt.figure()
plt.plot(x, y)
plt.show()

# f = open("results.txt", 'a')
# f.write(str(lr_asymmetry1))
# f.write('\t')
# f.write(str(lr_asymmetry2))
# f.write('\n')
# f.close()
