import matplotlib.pyplot as plt
import numpy as np

nk = 20
nb = 208

bands= np.loadtxt(f"work/bands.dat")

ans = np.reshape(bands, (nk,nb,4))

ans = np.array(ans)

for ib in range(nb):
  plt.plot(ans[:,ib,3], 'b')

plt.ylim((2,6.0))

plt.savefig("bands.png")
