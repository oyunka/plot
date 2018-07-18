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


fname = 'cep90/Ne.00.00.00.05.SPECTRUM/Ne.00.00.00.05_00000060.ang_spec'


pi = np.pi

with open(fname, "r") as fp:
    fp.readline()
    num_krad = int(fp.readline().strip("\n"))
    num_kang = int(fp.readline().strip("\n"))

data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None, skiprows=3)
wang_data = pd.read_csv("wang.txt", comment="#", delim_whitespace=True, header=None, skiprows=1)


krad = data.iloc[:,0]
kang = data.iloc[:,1]
kval = data.iloc[:,2]

cos_ang = wang_data.iloc[:,0]
kang_rad = wang_data.iloc[:,1]
weight = wang_data.iloc[:,2]

krad = krad.values.reshape([num_krad,num_kang])
kang = kang.values.reshape([num_krad, num_kang])
kval = kval.values.reshape([num_krad, num_kang])




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

"""
for i in range(0, num_krad_half):
    for j in range(0, num_kang):
        I2[j] += kval[i, j]
"""
I2[:] = kval[71, :]
plt.plot(kang_rad, I2,'.-')
plt.xlabel('kang_degree')

plt.show()
