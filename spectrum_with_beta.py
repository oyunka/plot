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

with open(fname, "r") as fp:
    fp.readline()
    num_krad = int(fp.readline().strip("\n"))
    num_kang = int(fp.readline().strip("\n"))

data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None, skiprows=3)
wang_data = pd.read_csv("wang.txt", comment="#", delim_whitespace=True, header=None, skiprows=1)
beta_data = pd.read_csv("beta_parameters.txt", comment="#", delim_whitespace=True, header=None)

print(beta_data)

print(wang_data)

krad = data.ix[:,0]
kang = data.ix[:,1]
kval = data.ix[:,2]

cos_ang = wang_data.ix[:,0]
kang_rad = wang_data.ix[:,1]
weight = wang_data.ix[:,2]


krad = np.reshape(krad, (num_krad, num_kang))
kang = np.reshape(kang, (num_krad, num_kang))
kval = np.reshape(kval, (num_krad, num_kang))


beta0 = beta_data.ix[:, 0]
beta1 = beta_data.ix[:, 1]
beta2 = beta_data.ix[:, 2]
beta3 = beta_data.ix[:, 3]
beta4 = beta_data.ix[:, 4]


def Legendre1():
    return  cos_ang

def Legendre2():
    first = np.zeros((num_kang), dtype = 'float')
    for m in range(0, num_kang):
        first[m] = 0.5 * (3*(cos_ang[m])**2 - 1)   
    return  first

def Legendre3():
    first = np.zeros((num_kang), dtype = 'float')

    for m in range(0, num_kang):
        first[m] =0.5 * (5 * pow(cos_ang[m], 3)  - 3 * cos_ang[m])
    return  first

def Legendre4():
    first = np.zeros((num_kang), dtype = 'float')

    for m in range(0, num_kang):
        first[m] = 1/8 * (35 * pow(cos_ang[m], 4) - 30 * pow (cos_ang[m], 2) + 3)
    return  first

legendre1 = Legendre1()
legendre2 = Legendre2()
legendre3 = Legendre3()
legendre4 = Legendre4()

print(legendre1)

I1 = np.empty((num_kang), dtype='float')
kang_degree = kang_rad * 57.2958
for i in range(0, num_kang):
    I1[i] = 1 + beta1 * legendre1[i] + beta2 * legendre2[i] + beta3 * legendre3[i] + beta4 * legendre4[i]

# plt.figure()
# plt.plot(kang_rad, I1)
# plt.show()

dk = np.zeros((num_krad), dtype = 'float')
weighted_kval = np.zeros((num_krad, num_kang), dtype = 'float')
momentum = krad[:,0]

for i in range(0,num_krad):
    if(i < num_krad):
        dk[i] = momentum[i] - momentum[i-1]
    if(i == 0) :
        dk[i] = momentum[i] - 0.0
I2 = np.zeros((num_kang), dtype='float')
num_krad_half = int(num_krad/2)


for i in range(0, num_krad_half):
    for j in range(0, num_kang):
        I2[j] += kval[i, j]*dk[i]* momentum[i] **2

# for i in range(0, num_krad_half):
#     for j in range(0, num_kang):
#         =  weighted_kval[i, j] *dk[i]* momentum[i] **2
#         print(dk[i]* momentum[i] **2 )

# print(I1,'\n',I2)

# plt.figure()
# plt.plot(kang_degree,I1, 'r--', kang_degree, I2, 'g^')
# plt.show()



# plt.subplot(2, 1, 1)
# plt.plot(kang_degree, I1, 'o-')
# plt.title('Plot with beta')

# plt.subplot(2, 1, 2)
# plt.plot(kang_degree, I2,'.-')
# plt.xlabel('kang_degree')

# plt.show()



fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_ylabel('Beta parameter plot', color=color)  # we already handled the x-label with ax1

plt.xlabel('Angle[°]')
ax1.plot(kang_degree, I2, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'

ax2.set_ylabel('Simulation result', color=color)  # we already handled the x-label with ax1
ax2.plot(kang_degree, I1, 'o')
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(fname+'_beta.pdf', bbox = 'tight')
plt.show()
