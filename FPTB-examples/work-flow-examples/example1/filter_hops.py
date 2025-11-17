import numpy as np

hops = np.loadtxt("hopping0.dat")

ans = []
for hop in hops:
  dx = int(hop[0])
  dy = int(hop[1])
  dz = int(hop[2])
  i = int(hop[3])
  j = int(hop[4])
  if (dx>=0 and (dx>0 or dy>0 or dz>0)) or ( j>=i and (dx==0 and dy==0 and dz==0)) :
    ans.append(hop)

ans = np.array(ans)

np.savetxt("hopping.dat",ans)
