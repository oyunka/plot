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
    cval = data.ix[:,2:].values
    print(cval.shape)
    #print(len(data.iloc[0,:]) - 2)

    krad = krad.reshape(num_krad, num_kang)
    kang = kang.reshape(num_krad, num_kang)

    #krad = np.reshape(krad, (num_krad, num_kang))
    #kang = np.reshape(kang, (num_krad, num_kang))
    #kval = np.reshape(kval, (num_krad, num_kang))

    norb = int(len(cval[0,:]) / 2)
    xang = np.cos(kang[0, :])

    lmax = 3
    mnum = 3
    mmax = 1
    val = np.zeros((norb, lmax+1, num_krad), dtype=np.complex)

    wdata = pd.read_csv(w_fname, comment="#", delim_whitespace=True, header=None, skiprows=0)
    #print(wdata)
    print("xang:")
    print(xang)
    wang = wdata.iloc[:, 2].values
    xang = wdata.iloc[:, 0].values

    print("# number of orbital : ", norb)
    m_list = [0, 0, -1, 0, 1]
    for iorb in range(norb):
        realp = cval[:,iorb*2].reshape(num_krad, num_kang)
        imagp = cval[:,iorb*2+1].reshape(num_krad, num_kang)
        #print(wang.shape)
        #print(imagp.shape)
        #quit()
        for irad in range(num_krad):
                for l in range(lmax+1):
                    m = m_list[iorb]
                    if np.abs(m) > l:
                        continue
                    assoc_leg = special.lpmv(np.abs(m), l, xang)
                    #fact = (-1)**( 0.5 * (m + np.abs(m))) * np.sqrt( (2*l+1) / (4*np.pi) * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                    """
                    #       ==============================
                    #       this part is missed in my code (td1c, surff_init)
                    """
                    fac = np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                    sph_val = fac * assoc_leg
                    #print(np.sum(sph_val  * sph_val * wang))
                    val[iorb, l, irad] = np.sum(wang * sph_val * (realp[irad,:] + IUNIT * imagp[irad,:]))

        print(type(val))
        print(val.shape)
    if 1:
        plt.figure()
        plt.title("first")

        pes_ang = np.zeros((num_krad, num_kang))

        for iorb in range(norb):
            realp = cval[:,iorb*2].reshape(num_krad, num_kang)
            imagp = cval[:,iorb*2+1].reshape(num_krad, num_kang)
            pes_ang += np.abs(realp + IUNIT * imagp)**2
            plt.plot(0.5*krad[:, 0]**2 * 27.2, krad[:, 0] * np.sum(np.abs(realp[:, :] + IUNIT*imagp[:, :])** 2 * wang, axis=1), "s-", label="orbital : {}, m : {}".format(iorb, m_list[iorb]))

        pes = np.zeros(num_krad)
        for irad in range(num_krad):
            pes[irad] += np.sum(pes_ang[irad,:] * wang)
        plt.plot(0.5*krad[:, 0]**2 * 27.2, pes * krad[:, 0], "+-", label = "total pes")

        #plt.xlim(10,25)
        plt.legend()
        plt.yscale("log")

    m_list2 = [-1, 0, 1]

    val_msum = np.zeros((mnum, lmax+1, num_krad), dtype=np.complex)
    amp = np.zeros((mnum, lmax+1, num_krad))
    phase = np.zeros((mnum, lmax+1, num_krad))
    for im, m in enumerate(m_list2):
        for iorb in range(norb):
            for l in range(lmax+1):
                if m == m_list[iorb]:
                    print("# m, iorb : ", m, ", ", iorb)
                    val_msum[im, l, :] += val[iorb, l, :]

    amp = np.abs(val_msum)
    #phase = np.angle(val_msum)
    #val_msum *= np.exp(IUNIT * np.pi * 0.5)
    #phase = np.arctan(np.imag(val_msum) / np.real(val_msum))
    #phase = np.real(val_msum) / np.abs(val_msum)

    #c = np.real(val_msum) / np.abs(val_msum)
    #s = np.imag(val_msum) / np.abs(val_msum)
    #phase = np.arctan2(s, c)
    phase = np.arctan2(np.imag(val_msum), np.real(val_msum))
    plt.show()

    plt.figure()
    plt.title("second")

    #ax1 = plt.subplot(211)
    #ax2 = plt.subplot(212)
    i = 0
    ax = []
    axt = []
    for im in range(mnum):
        for l in range(lmax+1):
            flag1 = (l == 0 and m_list2[im] ==  0)
            flag2 = (l == 1 and m_list2[im] == -1)
            flag3 = (l == 1 and m_list2[im] ==  0)
            flag4 = (l == 1 and m_list2[im] ==  1)
            flag5 = (l == 2 and m_list2[im] == -1)
            flag6 = (l == 2 and m_list2[im] ==  0)
            flag7 = (l == 2 and m_list2[im] ==  1)
            flag8 = (l == 3 and m_list2[im] ==  0)
            flag9 = (l == 3 and m_list2[im] == -1)
            flag10 = (l == 3 and m_list2[im] ==  1)
            #flag11 = (l == 2 and m_list2[im] ==  1)
            #flag12 = (l == 3 and m_list2[im] ==  0)
            #flag_tot = flag1 or
            flag_tot = flag1 or flag2  or flag3 or flag4 or flag5 or flag6 or flag7 or flag8 or flag9 or flag10

            if flag_tot:
                ax.append(plt.subplot(4,3,i+1))
                axt.append(ax[i].twinx())
                #print(i, len(ax))
                print(l, im)
                #v2 = val_msum / val_msum[2, 2, :]
                #ph = np.arctan(np.imag(v2) / np.real(v2))
                #ax[i].plot(0.5*krad[:, 0]**2 * 27.2, ph[im, l, :], "+", label = "l = {:+5d}, m = {:+5d}".format(l, m_list2[im]))

                ax[i].plot(0.5*krad[:, 0]**2 * 27.2, phase[im, l, :], "+-", label = "l = {:+5d}, m = {:+5d}".format(l, m_list2[im]))

                axt[i].plot(0.5*krad[:, 0]**2 * 27.2, amp[im, l, :], "k-", label = "l = {:+5d}, m = {:+5d}".format(l, m_list2[im]))
                #axt[i].set_yscale("log")
                #axt[i].set_ylim(10**-12, 0.035)

                ax[i].set_xlim(5, 15)
                #ax[i].set_xlim(10, 25)
                ax[i].set_ylim(-3.14, 3.14)

                ax[i].legend(loc='upper right')
                ax[i].grid()
                i = i + 1

    # plt.figure()
    #ax1 = plt.subplot(211)
    #ax2 = plt.subplot(212)
    i = 0
    ax = []
    axt = []
    for im in range(mnum):
        for l in range(lmax+1):
            flag1 = (l == 0 and m_list2[im] ==  0)
            flag2 = (l == 1 and m_list2[im] == -1)
            flag3 = (l == 1 and m_list2[im] ==  0)
            flag4 = (l == 1 and m_list2[im] ==  1)
            flag5 = (l == 2 and m_list2[im] == -1)
            flag6 = (l == 2 and m_list2[im] ==  0)
            flag7 = (l == 2 and m_list2[im] ==  1)
            flag8 = (l == 3 and m_list2[im] ==  0)
            flag9 = (l == 3 and m_list2[im] == -1)
            flag10 = (l == 3 and m_list2[im] ==  1)
            #flag11 = (l == 2 and m_list2[im] ==  1)
            #flag12 = (l == 3 and m_list2[im] ==  0)
            #flag_tot = flag1 or
            flag_tot = flag1 or flag2  or flag3 or flag4 or flag5 or flag6 or flag7 or flag8 or flag9 or flag10


            if flag_tot:
                ax.append(plt.subplot(4,3,i+1))
                axt.append(ax[i].twinx())
                #print(i, len(ax))
                print(l, im)
                v2 = val_msum / val_msum[2, 2, :]
                ph = np.arctan(np.imag(v2) / np.real(v2))
                ax[i].plot(0.5*krad[:, 0]**2 * 27.2, ph[im, l, :], "+", label = "l = {:+5d}, m = {:+5d}".format(l, m_list2[im]))

                #ax[i].plot(0.5*krad[:, 0]**2 * 27.2, phase[im, l, :], "+", label = "l = {:+5d}, m = {:+5d}".format(l, m_list2[im]))

                axt[i].plot(0.5*krad[:, 0]**2 * 27.2, amp[im, l, :], "k-", label = "l = {:+5d}, m = {:+5d}".format(l, m_list2[im]))
                #axt[i].set_yscale("log")
                #axt[i].set_ylim(10**-12, 0.035)


                ax[i].set_xlim(5, 15)
                #ax[i].set_xlim(10, 25)
                ax[i].set_ylim(-3.14, 3.14)
                ax[i].legend(loc='upper right')
                ax[i].grid()
                i = i + 1

    #plt.show()

    # calculation of norm
    norm = 0
    for im in range(len(m_list2)):
        for l in range(lmax):
            norm += np.sum(amp[im, l, :]**2)
    print(norm)

    norm = 0
    pes_ang = np.zeros((num_krad, num_kang))

    for iorb in range(norb):
        realp = cval[:,iorb*2].reshape(num_krad, num_kang)
        imagp = cval[:,iorb*2+1].reshape(num_krad, num_kang)
        pes_ang += np.abs(realp + IUNIT * imagp)**2

    pes = np.zeros(num_krad)
    for irad in range(num_krad):
        pes[irad] += np.sum(pes_ang[irad,:] * wang)

    print(np.sum(pes))


    # plotting
    targ_ene = 7.2
    ene = (0.5 * krad[:, 0]**2 * 27.2)
    targ_irad = np.argmin(np.abs(ene - targ_ene))
    print(targ_irad)
    #adpes = np.zeros((num_krad, ncos))

    nang = 24
    adpes0 = np.zeros(nang)
    adpes_sub = np.zeros((nang, len(m_list2)), dtype=np.complex)
    xang = np.linspace(0, np.pi, nang)
    xcos = np.cos(xang)
    print(xcos)
    #fig, ax = plt.subplots(5)


    for l in range(lmax+1):
        #if l == 1:
        #    continue
        for im, m in enumerate(m_list2):
            if np.abs(m) > l:
                continue
            else:
                assoc_leg = special.lpmv(np.abs(m), l, xcos)
                #fac = (-1)**( 0.5 * (m + np.abs(m)))  * np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                fac = np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                sph_val = fac * assoc_leg
                adpes_sub[:, im] += sph_val * val_msum[im, l, targ_irad]
    adpes0 += np.sum(np.abs(adpes_sub[:,:])**2, axis=1)*2

    adpes = np.zeros(num_kang, dtype=np.complex)
    adm = []
    adm.append(np.zeros(num_kang, dtype=np.complex))
    adm.append(np.zeros(num_kang, dtype=np.complex))
    adm.append(np.zeros(num_kang, dtype=np.complex))
    for im, m in enumerate(m_list2):
        for iorb in range(norb):
            if m == m_list[iorb]:
                realp = cval[:,iorb*2].reshape(num_krad, num_kang)[targ_irad, :]
                imagp = cval[:,iorb*2+1].reshape(num_krad, num_kang)[targ_irad, :]
                adm[im] += np.abs(realp + IUNIT * imagp)**2 * 2
                adpes += np.abs(realp + IUNIT * imagp)**2 * 2
    plt.figure()
    plt.plot(kang[0, :], adm[0], "+-", label="m = {}".format(m_list2[0]))
    plt.plot(kang[0, :], adm[1], "o-", label="m = {}".format(m_list2[1]))
    plt.plot(kang[0, :], adm[2], "x-", label="m = {}".format(m_list2[2]))
    plt.plot(kang[0, :], adpes, "s-",  label="total")
    plt.plot(xang, adpes0, "d-",  label="total from angular distribution")
    plt.plot(kang[0, :], adm[0]+adm[1]+adm[2], "d-", label="m total".format(m_list2[0]))

    plt.legend()
    plt.xlim(0, 3.14)



    adpes = np.zeros(nang)
    #fig, ax = plt.subplots(5)
    fig , ax = plt.subplots(4)
    for l in range(lmax+1):
        #if l == 1:
        #    continue
        for im, m in enumerate(m_list2):
            if np.abs(m) > l:
                continue
            else:
                assoc_leg = special.lpmv(np.abs(m), l, xcos)
                #fac = (-1)**( 0.5 * (m + np.abs(m)))  * np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                fac = np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
                sph_val = fac * assoc_leg
                adpes_sub[:, im] += sph_val * val_msum[im, l, targ_irad]
                ax[im].plot(xang, sph_val * np.abs(val_msum[im, l, targ_irad]), marker = (l+1, 1, 1), label="l = {}, m = {} ".format(l, m_list2[im]))
                ax[im].legend()
                ax[im].set_xlim(0, 3.14)
    adpes += np.sum(np.abs(adpes_sub[:,:])**2, axis=1)*2

    ax[3].plot(xang, adpes, "ko-", label="total")
    plt.legend()
    plt.xlim(0, 3.14)
    #plt.show()

    """
    adpes_orb = np.zeros((nang, norb), dtype=np.complex)
    for iorb in range(norb):
        for l in range(lmax+1):
            m = m_list[iorb]
            if np.abs(m) > l:
                continue
            assoc_leg = special.lpmv(np.abs(m), l, xcos)
            #fac = (-1)**( 0.5 * (m + np.abs(m))) * np.sqrt( (2*l+1) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
            fac = np.sqrt( (2.0*l+1.0) / 2.0 * np.math.factorial(l -  np.abs(m)) / np.math.factorial(l +  np.abs(m)))
            sph_val = fac * assoc_leg
            adpes_orb[:, iorb] += val[iorb, l, targ_irad] * sph_val

    adpes = np.sum(np.abs(adpes_orb[:,:])**2, axis=1)*2
    plt.figure()
    plt.plot(xang, adpes)
    plt.xlim(0, 3.14)
    """
    plt.show()
    return


if __name__ == "__main__":
    main()
