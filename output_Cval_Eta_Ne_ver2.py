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
Ev_const = 27.21138505

file_output = '19.1eV_Ne.0.0.5eV'

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
        print("one input is required.")
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
    cval = data.iloc[:,2:].values
    #print(cval.shape)
    #print(len(data.iloc[0,:]) - 2)

    krad = krad.reshape(num_krad, num_kang)
    kang = kang.reshape(num_krad, num_kang)


    norb = int(len(cval[0,:]) / 2)
    print('norb is', norb)
    xang = np.cos(kang[0, :])

    lmax = 3
    mnum = 3

    mmax = 1
    val = np.zeros((norb, lmax+1, num_krad), dtype=np.complex)

    wdata = pd.read_csv(w_fname, comment="#", delim_whitespace=True, header=None, skiprows=0)
    wang = wdata.iloc[:, 2].values
    xang = wdata.iloc[:, 0].values

    print("# number of orbital : ", norb)
    m_list = [0, 0, -1, 0, 1]
    for iorb in range(norb):
        realp = cval[:,iorb*2].reshape(num_krad, num_kang)
        imagp = cval[:,iorb*2+1].reshape(num_krad, num_kang)

        for irad in range(num_krad):
                for l in range(lmax+1):
                    m = m_list[iorb]
                    if np.abs(m) > l:
                        continue
                    assoc_leg = special.lpmv(np.abs(m), l, xang)
                    """
                    #       ==============================
                    #       this part is missed in my code (td1c, surff_init)
                    """
                    fac = np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                    sph_val = fac * assoc_leg
                    val[iorb, l, irad] = np.sum(wang * sph_val * (realp[irad,:] + IUNIT * imagp[irad,:]))
    m_list2 = [-1, 0, 1]
    val_msum = np.zeros((mnum, lmax+1, num_krad), dtype=np.complex)
    amp = np.zeros((mnum, lmax+1, num_krad))
    phase = np.zeros((mnum, lmax+1, num_krad))
    for im, m in enumerate(m_list2):
        for iorb in range(norb):
            for l in range(lmax+1):
                if m == m_list[iorb]:
#                     print("# im, m, iorb : ",im,", ",  m, ", ", iorb)
                    val_msum[im, l, :] += val[iorb, l, :]

    
    
    amp = np.abs(val_msum)
    phase = np.arctan2(np.imag(val_msum), np.real(val_msum))
    
#---------------------------------------------------------------------------------------------------
    Amp_sheet = np.zeros((mnum, lmax+1))
    Phase_sheet = np.zeros((mnum, lmax+1))

    for l in range(lmax+1): 
        for im in range(mnum):
            if(l == 0 and (im == 0 or im == 2 )):
                continue
           
            base = im
            v2 = val_msum / val_msum[base, 2, :]
            ph = np.arctan2(np.imag(v2), np.real(v2))
            
            col_amp = amp[im, l, 0:200]
            max_amp_pos = np.argmax(col_amp)
            col_phase = ph[im, l, 0:200]
            
            Amp_sheet[im, l] = np.amax(col_amp)
            Phase_sheet[im, l] = col_phase[max_amp_pos]
            
    col_name = ['l = 0', 'l = 1', 'l = 2', 'l = 3']
    df = pd.DataFrame(Amp_sheet, columns = col_name)
    df1 = pd.DataFrame(Phase_sheet, columns = col_name)
    df.index = df1.index =['m = -1', 'm = 0', 'm = 1']

    # df['l = 0'][0] = 'NaN'
    # df['l = 0'][2] = 'NaN'
    # df1['l = 0'][0] = 'NaN'
    # df1['l = 0'][2] = 'NaN'

    print(df)
    print(df1)
    
    df.to_csv(file_output + '_amp.txt', mode = 'w')
    df1.to_csv(file_output + '_phase.txt', mode = 'w')

    df_all = pd.concat([df, df1], axis = 1)

#             col.append(col_name)
#             df.index = 0.5*krad[60:80, 0]**2 * Ev_const
#             df.to_csv('amp_phase.csv', mode='w')
    df_all.to_csv('Amp_and_phase.csv', mode = 'w')

if __name__ == "__main__":
    main()
