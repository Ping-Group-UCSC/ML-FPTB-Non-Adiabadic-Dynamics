import matplotlib.pyplot as plt
import numpy as np

nk = 20
nb = 208
ans = []

for i in range(nk):
  evals = np.loadtxt(f"work/{i}/evals.dat")
  ans.append(evals)

ans = np.array(ans)

for ib in range(nb):
  plt.plot(ans[:,ib], 'b')


plt.ylim((2,6.0))

plt.savefig("bands.png")
