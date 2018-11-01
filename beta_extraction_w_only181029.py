import pandas as pd
import numpy as np
import sys
import math
import csv
import matplotlib.pyplot as plt
from scipy.special import*


# file_output = '14.3eV_Ne.0.0.5_w_only'
# file_output = '15.9eV_Ne.0.0.5_w_only'
#file_output = '19.1eV_Ne.0.0.5_w_only'
file_output = "tdhf"


argv = sys.argv
argc = len(argv)

if(argc == 2):
    fname = argv[1]
else:
    print("one input is required")

wang_data = pd.read_csv("wang.txt", comment="#", delim_whitespace=True, header=None, skiprows=1)
cos_ang = wang_data.iloc[:,0]
kang_rad = wang_data.iloc[:,1]
weight = wang_data.iloc[:,2]
def Beta(fname):
    data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None, skiprows=3)
    with open(fname, "r") as fp:
        fp.readline()
        num_krad = int(fp.readline().strip("\n"))
        num_kang = int(fp.readline().strip("\n"))

    krad = data.iloc[:,0]
    kang = data.iloc[:,1]
    kval = data.iloc[:,2]
    
    krad = krad.values.reshape([num_krad,num_kang])
    kang = kang.values.reshape([num_krad, num_kang])
    kval = kval.values.reshape([num_krad, num_kang])
    
    orbital_2s = 0
    momentum = krad[:,0]

    num_krad_half = int(num_krad/2)
    num_kang_half = int(num_kang/2)
    dk = np.zeros((num_krad), dtype = 'float')
    for i in range(0,num_krad):
        if(i < num_krad):
            dk[i] = momentum[i] - momentum[i-1]   
        if(i == 0) :
            dk[i] = momentum[i] - 0.0

    weighted_kval = np.zeros((num_krad), dtype = 'float')

    for i in range(0, num_krad):
        for j in range(0, num_kang):
            weighted_kval[i] += kval[i, j] * weight[j] * momentum[i] * momentum[i] * dk[i]

    i = np.argmax(weighted_kval)


    def Legendre0():
   
        all = 0
        for j in range(0, num_kang):
            all +=   kval[i, j]  * weight[j] * eval_legendre(0, cos_ang[j])

        return  all / 2
    def Legendre1():
        all = 0

        for j in range(0, num_kang):
            all +=   kval[i, j] * weight[j] * eval_legendre(1, cos_ang[j])

        return  all * 3/2

    def Legendre2():

        all = 0

        for j in range(0, num_kang):
            all +=   kval[i, j]  * weight[j] * eval_legendre(2, cos_ang[j])

        return  all * 5/2

    def Legendre3():
        all = 0

        for j in range(0, num_kang):
            all +=   kval[i, j] * weight[j] * eval_legendre(3, cos_ang[j])

        return all * 7/2

    def Legendre4():

        all = 0

        for j in range(0, num_kang):
            all +=   kval[i, j]  * weight[j] * eval_legendre(4, cos_ang[j])

        return all * 9/2

    beta0 = Legendre0()
    beta1 = Legendre1()/beta0
    beta2 = Legendre2()/beta0
    beta3 = Legendre3()/beta0
    beta4 = Legendre4()/beta0
    
    beta = pd.DataFrame([[beta1, beta2, beta3, beta4]], columns = ['beta1','beta2', 'beta3', 'beta4'])
    return   beta
betas = []

df = pd.DataFrame(index = [], columns = ['beta1','beta2', 'beta3', 'beta4'])


row = Beta(fname)
df = df.append(row)

df.index = ['cep0']

# print(df)

df1 = pd.concat([df, df], axis = 0)
print(df)

df.to_csv(file_output + '.csv', mode='w')
