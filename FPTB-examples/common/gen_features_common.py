import numpy as np
from scipy import special

class lookuptable:
  def __init__(self):
    self.dd = dict([])

  def set(self, lab1, u, lab2, v):
    #makes a mapping between lab1=u and lab2=v 
    if (lab1,lab2) not in self.dd:
      self.dd[(lab1,lab2)] = dict([])
      if lab1 != lab2:
        self.dd[(lab2,lab1)] = dict([])

    if u not in self.dd[(lab1,lab2)]:
      self.dd[(lab1,lab2)][u] = [v]
    else:
      self.dd[(lab1,lab2)][u].append(v)

    if v not in self.dd[(lab2,lab1)]:
      self.dd[(lab2,lab1)][v] = [u]
    else:
      self.dd[(lab2,lab1)][v].append(u)

  
  def get(self, *args):
    #args =  lab1, u1, lab2, u2 , ... , labn 
    # returns a list of v such that there are mappings between v and u1, v and u2, ...
    labn = args[-1]
    lab1 = args[0]
    u1 = args[1]
    if u1 not in self.dd[(lab1,labn)]:
      return []
    ans = self.dd[(lab1,labn)][u1]

    for i in range(2,len(args)-1,2):
      labm = args[i]
      um = args[i+1]
      if um not in self.dd[(labm,labn)]:
        return []
      ans = list(set(ans) & set(self.dd[(labm,labn)][um]))

    return ans



#Ewald summation
#cell[:,i] is the ith lattice vector.
#qs = charges
#ts[:,i] = positions of charge i
#sig = erf parameter
#rcut, gcut = realspace/recipspace cutoffs
def ewald_tote(cell,qs,ts,sig,rcut,gcut):
  shortrange = 0.0
  longrange = 0.0
  selfint = 0.0

  gvecs = 2*np.pi*np.linalg.inv(cell)
  vol = abs(np.linalg.det(cell))
  nq = len(qs)

  for ix in range(-rcut,rcut+1):
    for iy in range(-rcut,rcut+1):
      for iz in range(-rcut,rcut+1):
        for i in range(nq):
          for j in range(nq):

            if(i==j and ix==0 and iy==0 and iz==0):
              continue

            rdiff = np.linalg.norm(ts[:,i]-ts[:,j]+ np.matmul(cell,np.array([ix,iy,iz])))
            shortrange = shortrange + 0.5*qs[i]*qs[j]/rdiff * special.erfc(rdiff/np.sqrt(2.0)/sig)


  for gx in range(-gcut,gcut+1):
    for gy in range(-gcut,gcut+1):
      for gz in range(-gcut,gcut+1):

        if(gx==0 and gy==0 and gz==0):
          continue

        gvec = np.matmul(np.array([gx,gy,gz]),gvecs)

        strfac = 0.0

        for i in range(nq):
          strfac = strfac+qs[i]*np.exp(1.0j*np.dot(gvec,ts[:,i]))

        longrange = longrange + 2.0*np.pi/vol * np.exp(-0.5*sig**2 * np.linalg.norm(gvec)**2)/np.linalg.norm(gvec)**2 * np.linalg.norm(strfac)**2


  for i in range(nq):
    selfint = selfint + 1.0/np.sqrt(2.0*np.pi)/sig * qs[i]**2

  return shortrange+longrange-selfint

#ewald summation of potential (not total energy)
# i = index of atom to calculate potential at
# see micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf
def ewald_phi(i, cell,qs,ts,sig,rcut,gcut):
  shortrange = 0.0
  longrange = 0.0
  selfint = 0.0

  gvecs = 2*np.pi*np.linalg.inv(cell)
  vol = abs(np.linalg.det(cell))
  nq = len(qs)

  for ix in range(-rcut,rcut+1):
    for iy in range(-rcut,rcut+1):
      for iz in range(-rcut,rcut+1):
        for j in range(nq):

          if(i==j and ix==0 and iy==0 and iz==0):
            continue

          rdiff = np.linalg.norm(ts[:,i]-ts[:,j]+ np.matmul(cell,np.array([ix,iy,iz])))
          shortrange = shortrange + qs[j]/rdiff * special.erfc(rdiff/np.sqrt(2.0)/sig)

  for gx in range(-gcut,gcut+1):
    for gy in range(-gcut,gcut+1):
      for gz in range(-gcut,gcut+1):

        if(gx==0 and gy==0 and gz==0):
          continue

        gvec = np.matmul(np.array([gx,gy,gz]),gvecs)
        for j in range(nq):
          gr = np.dot(gvec,ts[:,i]-ts[:,j])
          longrange = longrange + 4.0*np.pi/vol * np.cos(gr) * np.exp(-0.5*sig**2 * np.linalg.norm(gvec)**2)/(np.linalg.norm(gvec)**2) * qs[j]


  selfint =  np.sqrt(2.0/np.pi)/sig * qs[i]

  return shortrange+longrange-selfint

