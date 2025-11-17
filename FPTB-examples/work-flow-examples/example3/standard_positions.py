import numpy as np

def posvec0(lut,atom1):
  #returns standard coordinates of an atom as a vector
#  ix iy iz for Cs
#  ix iy+0.5 iz+0.5 for Br1, etc
#  ix+0.5 iy+0.5 iz+0.5 for Pb
  ix = lut.get('atom',atom1,'ix')[0]
  iy = lut.get('atom',atom1,'iy')[0]
  iz = lut.get('atom',atom1,'iz')[0]
  typ = lut.get('atom',atom1,'atom_typ')[0]
  if typ=="Pb":
    t = np.array([0.5,0.5,0.5])
  elif typ=="Br1":
    t = np.array([0,0.5,0.5])
  elif typ=="Br2":
    t = np.array([0.5,0,0.5])
  elif typ=="Br3":
    t = np.array([0.5,0.5,0])
  elif typ=="Cs":
    t = np.array([0,0,0])
  return np.array([ix,iy,iz])+t

#translate atoms to be centered on standard positions:
# LX0,LY0,LZ0 = unitcell dimensions 
def translate_posvec0(lut, atom, c_old, LX0, LY0, LZ0):
  uc = np.array([LX0,LY0,LZ0])
  c_new = c_old - uc*posvec0(lut,atom) + 0.5*uc
  c_new = np.mod(c_new, uc)
  c_new = c_new + uc*posvec0(lut,atom) - 0.5*uc
  return c_new


