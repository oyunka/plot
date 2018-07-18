import numpy as np
import pandas as pd
import sys
import csv
import matplotlib.pyplot as plt

fname = sys.argv[1]
df = pd.read_csv(sys.argv[1], index_col = None, header = None, delim_whitespace = True, comment = "#")
data = df.as_matrix()[:, :]
time = data[:, 1]
efield = data[:, 3]


newt = np.linspace(time[0], (time[1] - time[0]) * len(time) * 11,  len(time) * 11)
efield = np.concatenate([efield, np.zeros(10 * len(time))])
freq = np.fft.fftfreq(len(newt), time[1] - time[0])
fft = np.fft.fft(efield)


df2 = pd.DataFrame(np.array([np.abs(freq), np.abs(fft), np.real(fft), np.imag(fft), (np.abs(fft))**2]).T)

df2.to_csv("fft.txt", index=None, header=None,
           sep="\t",
           float_format="%20.10E",
           quoting=csv.QUOTE_NONE,
           #quotechar=' ',
           #quotechar='"',
           #quotechar="'",
           #escapechar=" ",
           escapechar=""
           )
<<<<<<< HEAD
=======

print("what is going on")
print("calc-divide")
>>>>>>> calc-divide
