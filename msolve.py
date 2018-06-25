import pandas as pd
import numpy as np
from matplotlib import pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import sys

from scipy.interpolate import lagrange
from scipy.interpolate import interp1d
from scipy.interpolate import spline
from scipy.interpolate import BSpline
from scipy import special
argv = sys.argv
argc = len(argv)
IUNIT = 0 + 1j

def main():

    w_fname = "wang.txt"
    
    if(argc == 2):
        fname = argv[1]#"H_ecs_r40toinf_th15_s2t20.surff"
        #outname = fname
        outname = argv[1].split("/")[-1].strip(".pdf")
    elif(argc == 3):
        fname = argv[1]
        outname = argv[2].strip(".pdf")
    else:
        printf("one input is required.")
        quit()


        
    with open(fname, "r") as fp:
        line = fp.readline()
        if len(list(line.strip()))  != 0:
            if list(line.strip())[0] == "#":
                #print(line)
                num_krad = int(fp.readline().strip("\n"))
                num_kang = int(fp.readline().strip("\n"))
                skipnum = 3
            else:            
                #print("a")
                skipnum = 2
                num_krad = int(line.strip("\n"))
                num_kang = int(fp.readline().strip("\n"))

    data = pd.read_csv(fname, comment="#", delim_whitespace=True, header=None, skiprows=skipnum)
    krad = data.iloc[:,0].values
    kang = data.iloc[:,1].values
    cval = data.ix[:,2:].values    

    
    krad = krad.reshape(num_krad, num_kang)
    kang = kang.reshape(num_krad, num_kang)
    

    norb = int(len(cval[0,:]) / 2)
    xang = np.cos(kang[0, :])

    m_list = [0, 0, -1, 0, 1]
    m_list2 = [-1, 0, 1]
    targ_ene = 7.2
    ene = (0.5 * krad[:, 0]**2 * 27.2)
    targ_irad = np.argmin(np.abs(ene - targ_ene))
    print(targ_irad)    
    

    adpes = np.zeros(num_kang)
    adm = []
    adm.append(np.zeros(num_kang))
    adm.append(np.zeros(num_kang))
    adm.append(np.zeros(num_kang))
    for im, m in enumerate(m_list2):
        for iorb in range(norb):
            if m == m_list[iorb]:
                realp = cval[:,iorb*2].reshape(num_krad, num_kang)[targ_irad, :]
                imagp = cval[:,iorb*2+1].reshape(num_krad, num_kang)[targ_irad, :]
                adm[im] += np.abs(realp + IUNIT * imagp)**2 * 2
                adpes += np.abs(realp + IUNIT * imagp)**2 * 2

    if 0:
        plt.figure()
        plt.plot(kang[0, :], adm[0], "+-", label="m = {}".format(m_list2[0]))
        plt.plot(kang[0, :], adm[1], "o-", label="m = {}".format(m_list2[1]))
        plt.plot(kang[0, :], adm[2], "x-", label="m = {}".format(m_list2[2]))
        plt.plot(kang[0, :], adpes, "s-",  label="total")
        plt.plot(kang[0, :], adm[0]+adm[1]+adm[2], "d-", label="m total".format(m_list2[0]))
        
        plt.legend()
        plt.xlim(0, 3.14)
        plt.show()

    with open("data.txt", "w") as f:
        f.write("# {} {} {} \n".format(*m_list2))
        for iang in range(num_kang):
            for i, m in enumerate(m_list2):
                f.write(" {: 20.10E} ".format(adm[i][iang]))
            f.write("\n")
    return


if __name__ == "__main__":
    main()
