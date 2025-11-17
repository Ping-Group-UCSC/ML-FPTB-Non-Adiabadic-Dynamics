import numpy as np
import os
import gen_features_common as gfc


rcs=2
rc=10
# count up to m-th neighbor shells
mshells_pbpb = 2
mshells_pbi = 2
mshells_ipb = 2
mshells_ii = 2
mshells_pbma = 1
mshells_ima = 1
mshells0 = 2 # this is just for constructing some hoppings like nnii_cage, not for making features.
mshells_max = max(mshells_pbpb,mshells_pbi,mshells_ipb,mshells_ii,mshells_pbma,mshells_ima,mshells0)

# Load atom coordinates
#

NN = 96 # 2x2x2 supercell, all atoms
NX = 2
NY = 2
NZ = 2
LX = 12.532 #supercell dimensions
LY = 12.657
LZ = 12.532
LX0 = LX/NX #unitcell dimensions
LY0 = LY/NY
LZ0 = LZ/NZ


## max/min 
pb_pb_min=110
pb_i_min=110
i_i_min=110
i_pb_min=110 
pb_pb_max=0
i_i_max=0
pb_i_max=0
i_pb_max=0

LL=np.array([LX,LY,LZ])
LL=np.reshape(LL, [-1, 1, 3])
epsilon = 1e-14
f_coord = open('./coord.xyz')
coord = []
#sym=[]
while f_coord.readline():
  f_coord.readline()
  for ii in range(NN):
    line_in = f_coord.readline()
#    sym.append(line_in.split()[0])
    coord.append(np.reshape(np.array(line_in.split()[1:]).astype(np.float),[-1,1]))
coord = np.reshape(np.concatenate(coord, axis = 0), [-1, NN, 1, 3])



#coord = np.zeros((5,NN,1,3))


import sys 
nframe = np.shape(coord)[0]
print ("Total number of frames %d\n"%(nframe))
sys.stdout.flush()

# > r_ij, a (nframe, NN, NN) tensor, (n, i, j) means for frame n, the distance between i and j.
# 
# > coord_ij,a (nframe, NN, NN, 4) tensor, (n, i, j, k) means for frame n, the k-th element of sij, where sij=(xij,yij,zij,rij)/rij^2
# 
# calculate the rij matrix
nframe = np.shape(coord)[0]
coord_ij = np.zeros((nframe, NN, NN, 3))
for ii in range(NN):
  coord_ij[:,ii,:,:3] = np.mod(coord[:,:,0,:] - np.reshape(coord[:,ii,0,:], [-1, 1, 3]) + LL, LL) 


#for jj in range(NN):
#  out=open (f"new{jj}.xyz",'w')
#  out.write(f"{NN} \n\n")
#  for ii in range(NN):
#    out.write(f"{sym[ii]} {coord_ij[0,ii,jj,:][:]}\n")
#  out.close()
    
r_ij =  np.linalg.norm(coord_ij, axis=3, keepdims=True)
#coord_ij = np.concatenate([coord_ij, r_ij], axis = 3) / (r_ij+EPSILON)**2
coord_ij = np.concatenate([coord_ij, r_ij], axis = 3) #/ (r_ij+EPSILON)**2


#Sort j in coord_ij according to the magnitude of rij
sort_idx = np.argsort(r_ij, axis=2)
np.shape(sort_idx)
r_ij_sorted = np.zeros(np.shape(r_ij))
coord_ij_sorted = np.zeros(np.shape(coord_ij))
for nn in range(nframe):
  for ii in range(NN):
    for jj in range(NN):
      r_ij_sorted[nn,ii,jj] = r_ij[nn,ii,sort_idx[nn,ii,jj]]
      coord_ij_sorted[nn,ii,jj,:] = coord_ij[nn,ii,sort_idx[nn,ii,jj],:]


def smooth_rij(rij):
#  return (1.0/rij)
  if rij < rcs :
    return 1.0/rij
  elif rij > rc :
    return 0
  else:
    return ((1.0/rij)*(0.5*np.cos(np.pi*(rij-rcs)/(rc-rcs))+0.5))
      


#
# initialize lookuptables
#
def makelut(atomindexdat, orbindexdat):
  lut = lookuptable()
  atom_index_file = open(atomindexdat,'r')
  for line in atom_index_file:
    i,atom_typ,ix,iy,iz = line.split()
    ## WAS
    i=int(i)-1
    ###
    lut.set("atom", int(i), "atom_typ", atom_typ)
    lut.set("atom", int(i), "ix", int(ix))
    lut.set("atom", int(i), "iy", int(iy))
    lut.set("atom", int(i), "iz", int(iz))
  atom_index_file.close()

  orb_index_file = open(orbindexdat,'r')
  for line in orb_index_file:
    ll = line.split()
    orb_idx = int(ll[0])
    orb_typ = ll[1]
    ix = int(ll[2])
    iy = int(ll[3])
    iz = int(ll[4])
    atom_typ = orb_typ.split(':')[0]
  
    lut.set("orb", orb_idx, "orb_typ", orb_typ)
  
    if orb_typ=="Pb:s":
      spin = ll[5]
      lut.set("orb", orb_idx, "spin", spin)

    if orb_typ=="Pb:p":
      pdir = ll[5]
      spin = ll[6]
      lut.set("orb", orb_idx, "pdir", pdir)
      lut.set("orb", orb_idx, "spin", spin)

    if orb_typ=="I:p":
      ddir = ll[5]
      pdir = ll[6]
      spin = ll[7]
      lut.set("orb", orb_idx, "dir", ddir)
      lut.set("orb", orb_idx, "pdir", pdir)
      lut.set("orb", orb_idx, "spin", spin)
    
      if ddir=='x':
        atom_typ="I1"
      if ddir=='y':
        atom_typ="I2"
      if ddir=='z':
        atom_typ="I3"
    
    lut.set("orb", orb_idx, "atom_typ", atom_typ)
    lut.set("orb", orb_idx, "ix", ix)
    lut.set("orb", orb_idx, "iy", iy)
    lut.set("orb", orb_idx, "iz", iz)

    atom_idx = lut.get("atom_typ",atom_typ,"ix",ix,"iy",iy,"iz",iz,"atom")[0]
    lut.set("orb", orb_idx, "atom", atom_idx)


  orb_index_file.close()
  return lut




def posvec0(lut,atom1):
  #returns canonical coordinates of an atom as a vector
  ix = lut.get('atom',atom1,'ix')[0]
  iy = lut.get('atom',atom1,'iy')[0]
  iz = lut.get('atom',atom1,'iz')[0]
  typ = lut.get('atom',atom1,'atom_typ')[0]
  if typ=="Pb":
    t = np.array([0.0,0.0,0.0])
  elif typ=="I1":
    t = np.array([0.5,0,0])
  elif typ=="I2":
    t = np.array([0,0.5,0])
  elif typ=="I3":
    t = np.array([0,0,0.5])
  elif typ=="N":
    t = np.array([0.5,0.5,0.5])
  return np.array([ix,iy,iz])+t


lut = gfc.makelut("atom_index.dat","orb_index.dat")

#load shells of n-th neighbor atoms
pbpb_shells = np.loadtxt("shells/pbpb_shells.dat",dtype=int)
pbi_shells = np.loadtxt("shells/pbi_shells.dat",dtype=int)
pbma_shells = np.loadtxt("shells/pbma_shells.dat",dtype=int)
i1pb_shells = np.loadtxt("shells/i1pb_shells.dat",dtype=int)
i2pb_shells = np.loadtxt("shells/i2pb_shells.dat",dtype=int)
i3pb_shells = np.loadtxt("shells/i3pb_shells.dat",dtype=int)
i1ma_shells = np.loadtxt("shells/i1ma_shells.dat",dtype=int)
i2ma_shells = np.loadtxt("shells/i2ma_shells.dat",dtype=int)
i3ma_shells = np.loadtxt("shells/i3ma_shells.dat",dtype=int)
i1i_shells = np.loadtxt("shells/i1i_shells.dat",dtype=int)
i2i_shells = np.loadtxt("shells/i2i_shells.dat",dtype=int)
i3i_shells = np.loadtxt("shells/i3i_shells.dat",dtype=int)
ipb_shells = [i1pb_shells,i2pb_shells,i3pb_shells]
ima_shells = [i1ma_shells,i2ma_shells,i3ma_shells]
ii_shells = [i1i_shells,i2i_shells,i3i_shells]

n_pbpb_shells = np.loadtxt("shells/n_pbpb_shells.dat",dtype=int)
n_pbi_shells = np.loadtxt("shells/n_pbi_shells.dat",dtype=int)
n_pbma_shells = np.loadtxt("shells/n_pbma_shells.dat",dtype=int)
n_i1pb_shells = np.loadtxt("shells/n_i1pb_shells.dat",dtype=int)
n_i2pb_shells = np.loadtxt("shells/n_i2pb_shells.dat",dtype=int)
n_i3pb_shells = np.loadtxt("shells/n_i3pb_shells.dat",dtype=int)
n_i1ma_shells = np.loadtxt("shells/n_i1ma_shells.dat",dtype=int)
n_i2ma_shells = np.loadtxt("shells/n_i2ma_shells.dat",dtype=int)
n_i3ma_shells = np.loadtxt("shells/n_i3ma_shells.dat",dtype=int)
n_i1i_shells = np.loadtxt("shells/n_i1i_shells.dat",dtype=int)
n_i2i_shells = np.loadtxt("shells/n_i2i_shells.dat",dtype=int)
n_i3i_shells = np.loadtxt("shells/n_i3i_shells.dat",dtype=int)
n_ipb_shells = [n_i1pb_shells,n_i2pb_shells,n_i3pb_shells]
n_ima_shells = [n_i1ma_shells,n_i2ma_shells,n_i3ma_shells]
n_ii_shells = [n_i1i_shells,n_i2i_shells,n_i3i_shells]


for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):

      #print ("%s %s %s" % (ix,iy,iz))
      atom1 = lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom" )[0]

      #pb-pb shells
      for ishell in range(mshells_max): 
        j1shell = n_pbpb_shells[ishell+1] 
        j2shell = n_pbpb_shells[ishell+2] 

        for jshell in range(j1shell,j2shell):
          dx,jx = np.divmod( ix + pbpb_shells[jshell,0] , NX)
          dy,jy = np.divmod( iy + pbpb_shells[jshell,1] , NY)
          dz,jz = np.divmod( iz + pbpb_shells[jshell,2] , NZ)
          atom2 = lut.get("atom_typ","Pb","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]

          lut.set("atom", atom1, "pbpb_shell_%s" % (ishell+1), (atom2,dx,dy,dz))

      #pb-i shells
      for ishell in range(mshells_max): 
        j1shell = n_pbi_shells[ishell+0] 
        j2shell = n_pbi_shells[ishell+1] 

        for jshell in range(j1shell,j2shell):
          dx,jx = np.divmod( ix + pbi_shells[jshell,0] , NX)
          dy,jy = np.divmod( iy + pbi_shells[jshell,1] , NY)
          dz,jz = np.divmod( iz + pbi_shells[jshell,2] , NZ)
          jtyp = pbi_shells[jshell,3]
          atom2 = lut.get("atom_typ","I%s" % jtyp,"ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]

          lut.set("atom", atom1, "pbi_shell_%s" % (ishell+1), (atom2,dx,dy,dz))

      #pb-ma shells
      for ishell in range(mshells_max): 
        j1shell = n_pbma_shells[ishell+0] 
        j2shell = n_pbma_shells[ishell+1] 

        for jshell in range(j1shell,j2shell):
          dx,jx = np.divmod( ix + pbma_shells[jshell,0] , NX)
          dy,jy = np.divmod( iy + pbma_shells[jshell,1] , NY)
          dz,jz = np.divmod( iz + pbma_shells[jshell,2] , NZ)
          atom2 = lut.get("atom_typ","N","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]

          lut.set("atom", atom1, "pbma_shell_%s" % (ishell+1), (atom2,dx,dy,dz))

      
      for ityp in range(1,4):
        atom1 = lut.get("atom_typ","I%s" % ityp,"ix",ix+1,"iy",iy+1,"iz",iz+1,"atom" )[0]
        #i-pb shells
        for ishell in range(mshells_max): 
          j1shell = n_ipb_shells[ityp-1][ishell+0] 
          j2shell = n_ipb_shells[ityp-1][ishell+1] 

          for jshell in range(j1shell,j2shell):
            dx,jx = np.divmod( ix + ipb_shells[ityp-1][jshell,0] , NX)
            dy,jy = np.divmod( iy + ipb_shells[ityp-1][jshell,1] , NY)
            dz,jz = np.divmod( iz + ipb_shells[ityp-1][jshell,2] , NZ)
            atom2 = lut.get("atom_typ","Pb","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]

            lut.set("atom", atom1, "ipb_shell_%s" % (ishell+1), (atom2,dx,dy,dz))

        #i-ma shells
        for ishell in range(mshells_max): 
          j1shell = n_ima_shells[ityp-1][ishell+0] 
          j2shell = n_ima_shells[ityp-1][ishell+1] 

          for jshell in range(j1shell,j2shell):
            dx,jx = np.divmod( ix + ima_shells[ityp-1][jshell,0] , NX)
            dy,jy = np.divmod( iy + ima_shells[ityp-1][jshell,1] , NY)
            dz,jz = np.divmod( iz + ima_shells[ityp-1][jshell,2] , NZ)
            atom2 = lut.get("atom_typ","N","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]

            lut.set("atom", atom1, "ima_shell_%s" % (ishell+1), (atom2,dx,dy,dz))

        #i-i shells
        for ishell in range(mshells_max): 
          j1shell = n_ii_shells[ityp-1][ishell+1] 
          j2shell = n_ii_shells[ityp-1][ishell+2] 

          for jshell in range(j1shell,j2shell):

            dx,jx = np.divmod( ix + ii_shells[ityp-1][jshell,0] , NX)
            dy,jy = np.divmod( iy + ii_shells[ityp-1][jshell,1] , NY)
            dz,jz = np.divmod( iz + ii_shells[ityp-1][jshell,2] , NZ)
            jtyp = ii_shells[ityp-1][jshell,3]
            atom2 = lut.get("atom_typ","I%s" % jtyp,"ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]

            lut.set("atom", atom1, "ii_shell_%s" % (ishell+1), (atom2,dx,dy,dz))
 

#translate atoms to be centered on canonical positions:
#  ix iy iz for Pb
#  ix+0.5 iy iz for I1, etc
#  ix+0.5 iy+0.5 iz+0.5 for MA
for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):

      atom = lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for nn in range(nframe):
        coord[nn,atom,0,:] = coord[nn,atom,0,:] + np.array([LX0/2,LY0/2,LZ0/2])
        coord[nn,atom,0,:] = np.mod(coord[nn,atom,0,:], np.array([LX0,LY0,LZ0])) 
        coord[nn,atom,0,:] = coord[nn,atom,0,:] - np.array([LX0/2,LY0/2,LZ0/2]) + np.array([ix*LX0,iy*LY0,iz*LZ0])

      atom = lut.get("atom_typ","I1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for nn in range(nframe):
        coord[nn,atom,0,:] = coord[nn,atom,0,:] + np.array([0,LY0/2,LZ0/2])
        coord[nn,atom,0,:] = np.mod(coord[nn,atom,0,:], np.array([LX0,LY0,LZ0])) 
        coord[nn,atom,0,:] = coord[nn,atom,0,:] - np.array([0,LY0/2,LZ0/2]) + np.array([ix*LX0,iy*LY0,iz*LZ0])

      atom = lut.get("atom_typ","I2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for nn in range(nframe):
        coord[nn,atom,0,:] = coord[nn,atom,0,:] + np.array([LX0/2,0,LZ0/2])
        coord[nn,atom,0,:] = np.mod(coord[nn,atom,0,:], np.array([LX0,LY0,LZ0])) 
        coord[nn,atom,0,:] = coord[nn,atom,0,:] - np.array([LX0/2,0,LZ0/2]) + np.array([ix*LX0,iy*LY0,iz*LZ0])

      atom = lut.get("atom_typ","I3","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for nn in range(nframe):
        coord[nn,atom,0,:] = coord[nn,atom,0,:] + np.array([LX0/2,LY0/2,0])
        coord[nn,atom,0,:] = np.mod(coord[nn,atom,0,:], np.array([LX0,LY0,LZ0])) 
        coord[nn,atom,0,:] = coord[nn,atom,0,:] - np.array([LX0/2,LY0/2,0]) + np.array([ix*LX0,iy*LY0,iz*LZ0])

      for atom_typ in ["N","C","H_N1","H_N2","H_N3","H_C1","H_C2","H_C3"]:
        atom = lut.get("atom_typ",atom_typ,"ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
        for nn in range(nframe):
          coord[nn,atom,0,:] = np.mod(coord[nn,atom,0,:], np.array([LX0,LY0,LZ0])) 
          coord[nn,atom,0,:] = coord[nn,atom,0,:] + np.array([ix*LX0,iy*LY0,iz*LZ0])




#
# construct hoppings of each type here
# dx, dy, dz, orb1, orb2

pbp = []
for i in lut.get("orb_typ","Pb:p","spin","up","orb"):
  pbp.append([0,0,0,i,i])

pbs = []
for i in lut.get("orb_typ","Pb:s","spin","up","orb"):
  pbs.append([0,0,0,i,i])

ip_para = []
ip_perp = []
for orb1 in lut.get("orb_typ","I:p","spin","up","orb"):
  ddir = lut.get("orb",orb1,"dir")[0]
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  if ddir==pdir:
    ip_para.append([0,0,0,orb1,orb1])
  else:
    ip_perp.append([0,0,0,orb1,orb1])


pppi = []
ppsig= []
spsig= []

ppn_i_to_pb = []  #used to be called jn_i_to_pb
ppn_pb_to_i = []  #used to be called jn_pb_to_i
spn_pb_to_i = []

for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):

      atom1 = lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for neighbor in lut.get("atom", atom1, "pbi_shell_1"): 
        atom2,dx,dy,dz = neighbor

        for orb1 in lut.get("atom",atom1,"spin","up","orb"):
          for orb2 in lut.get("atom",atom2,"spin","up","orb"):

            ddir2 = lut.get("orb",orb2,"dir")[0]
            pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
            if lut.get("orb",orb1,"orb_typ")[0] == "Pb:s":
              if ddir2==pdir2:
                spsig.append([dx,dy,dz,orb1,orb2])
              else:
                spn_pb_to_i.append([dx,dy,dz,orb1,orb2])
            elif lut.get("orb",orb1,"orb_typ")[0] == "Pb:p":
              pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
              if ddir2==pdir1:
                if ddir2==pdir2:
                  ppsig.append([dx,dy,dz,orb1,orb2])
                else:
                  ppn_pb_to_i.append([dx,dy,dz,orb1,orb2])
              else:
                if pdir1==pdir2:
                  pppi.append([dx,dy,dz,orb1,orb2])
                elif ddir2==pdir2:
                  ## WAS 
#                  ppn_i_to_pb.append([dx,dy,dz,orb1,orb2])
                  ## WAS to be consisent with old data 
                  ppn_i_to_pb.append([-dx,-dy,-dz,orb2,orb1])
                  ###


nnpbpb = [] # sigma bonds involving p orbitals only
for ix in range(NX): 
  for iy in range(NY):
    for iz in range(NZ):

      atom1 = lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      
      for neighbor in lut.get("atom", atom1, "pbpb_shell_1"): 
        atom2,dx,dy,dz = neighbor

        jx = lut.get("atom",atom2,"ix")[0]-1
        jy = lut.get("atom",atom2,"iy")[0]-1
        jz = lut.get("atom",atom2,"iz")[0]-1
        
        if jx+NX*dx >= ix and jy+NY*dy >= iy and jz+NZ*dz >= iz:

          if jx+NX*dx > ix:
            bdir="x"
          elif jy+NY*dy > iy:
            bdir="y"
          elif jz+NZ*dz > iz:
            bdir="z"

          for orb1 in lut.get("atom",atom1,"spin","up","orb_typ","Pb:p","orb"):
            for orb2 in lut.get("atom",atom2,"spin","up","orb_typ","Pb:p","orb"):
              pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
              pdir2 = lut.get("orb",orb2,"pdir")[0].lower()

              if pdir1==pdir2 and bdir==pdir1 :
                nnpbpb.append([dx,dy,dz,orb1,orb2])


nnii_cage = [] #includes ii nearest neighbor bonds involving p-orbitals perpendicular to I-Pb 
for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):

      atom1 = lut.get("atom_typ","I1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for neighbor in lut.get("atom", atom1, "ii_shell_1"): 
        atom2,dx,dy,dz = neighbor

        for orb1 in lut.get("atom",atom1,"spin","up","orb"):
          for orb2 in lut.get("atom",atom2,"spin","up","orb"):
            pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
            pdir2 = lut.get("orb",orb2,"pdir")[0].lower()

            ddir1 = lut.get("orb",orb1,"dir")[0]
            ddir2 = lut.get("orb",orb2,"dir")[0]
            if pdir1!=pdir2 and pdir1==ddir2 and pdir2==ddir1:
              nnii_cage.append([dx,dy,dz,orb1,orb2])

      atom1 = lut.get("atom_typ","I2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      for neighbor in lut.get("atom", atom1, "ii_shell_1"): 
        atom2,dx,dy,dz = neighbor
        atom_typ2 = lut.get("atom",atom2,"atom_typ")[0]
        if atom_typ2 == "I3":

          for orb1 in lut.get("atom",atom1,"spin","up","orb"):
            for orb2 in lut.get("atom",atom2,"spin","up","orb"):
              pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
              pdir2 = lut.get("orb",orb2,"pdir")[0].lower()

              ddir1 = lut.get("orb",orb1,"dir")[0]
              ddir2 = lut.get("orb",orb2,"dir")[0]
              if pdir1!=pdir2 and pdir1==ddir2 and pdir2==ddir1:
                nnii_cage.append([dx,dy,dz,orb1,orb2])


nnii_sig = [] #sigma bonding "across" I-Pb-I only
for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):

      atom1 = lut.get("atom_typ","I1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      atom2 = lut.get("atom_typ","I2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
      atom3 = lut.get("atom_typ","I3","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]

      for atom4 in [atom1,atom2,atom3]:

        for neighbor in lut.get("atom", atom4, "ii_shell_2"): 
          atom5,dx,dy,dz = neighbor

          jx = lut.get("atom",atom5,"ix")[0]-1
          jy = lut.get("atom",atom5,"iy")[0]-1
          jz = lut.get("atom",atom5,"iz")[0]-1

          if jx+NX*dx >= ix and jy+NY*dy >= iy and jz+NZ*dz >= iz:

            if jx+NX*dx > ix:
              bdir="x"
            elif jy+NY*dy > iy:
              bdir="y"
            elif jz+NZ*dz > iz:
              bdir="z"

            for orb4 in lut.get("atom",atom4,"spin","up","orb_typ","I:p","orb"):
              for orb5 in lut.get("atom",atom5,"spin","up","orb_typ","I:p","orb"):
                pdir4 = lut.get("orb",orb4,"pdir")[0].lower()
                pdir5 = lut.get("orb",orb5,"pdir")[0].lower()

                ddir4 = lut.get("orb",orb4,"dir")[0]
                ddir5 = lut.get("orb",orb5,"dir")[0]
                if pdir4==pdir5 and pdir4==bdir and pdir4==ddir4 :
                  nnii_sig.append([dx,dy,dz,orb4,orb5])


socpb = [] 

soci_par = [] #this is when one of the orbitals is parallel to the bond 
soci_perp = [] #this is when both of the orbitals are perpendicular to the bond

for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):

      for atom1 in lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):

        orb1 = lut.get("atom",atom1,"spin","up","orb_typ","Pb:p","pdir","X","orb")[0]
        orb2 = lut.get("atom",atom1,"spin","dn","orb_typ","Pb:p","pdir","X","orb")[0]
        orb3 = lut.get("atom",atom1,"spin","up","orb_typ","Pb:p","pdir","Y","orb")[0]
        orb4 = lut.get("atom",atom1,"spin","dn","orb_typ","Pb:p","pdir","Y","orb")[0]
        orb5 = lut.get("atom",atom1,"spin","up","orb_typ","Pb:p","pdir","Z","orb")[0]
        orb6 = lut.get("atom",atom1,"spin","dn","orb_typ","Pb:p","pdir","Z","orb")[0]

        socpb.append([0,0,0,orb1,orb3])
        socpb.append([0,0,0,orb2,orb4])
        socpb.append([0,0,0,orb2,orb5])
        socpb.append([0,0,0,orb1,orb6])
        socpb.append([0,0,0,orb3,orb6])
        socpb.append([0,0,0,orb4,orb5])

      for atom2 in lut.get("atom_typ","I1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):

        orb1 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","X","orb")[0]
        orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","X","orb")[0]
        orb3 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","Y","orb")[0]
        orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","Y","orb")[0]
        orb5 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","Z","orb")[0]
        orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","Z","orb")[0]

        soci_par.append([0,0,0,orb1,orb3])
        soci_par.append([0,0,0,orb2,orb4])
        soci_par.append([0,0,0,orb2,orb5])
        soci_par.append([0,0,0,orb1,orb6])
        soci_perp.append([0,0,0,orb3,orb6])
        soci_perp.append([0,0,0,orb4,orb5])

      for atom2 in lut.get("atom_typ","I2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):

        orb1 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","X","orb")[0]
        orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","X","orb")[0]
        orb3 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","Y","orb")[0]
        orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","Y","orb")[0]
        orb5 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","Z","orb")[0]
        orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","Z","orb")[0]

        soci_par.append([0,0,0,orb1,orb3])
        soci_par.append([0,0,0,orb2,orb4])
        soci_perp.append([0,0,0,orb2,orb5])
        soci_perp.append([0,0,0,orb1,orb6])
        soci_par.append([0,0,0,orb3,orb6])
        soci_par.append([0,0,0,orb4,orb5])

      for atom2 in lut.get("atom_typ","I3","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):

        orb1 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","X","orb")[0]
        orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","X","orb")[0]
        orb3 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","Y","orb")[0]
        orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","Y","orb")[0]
        orb5 = lut.get("atom",atom2,"spin","up","orb_typ","I:p","pdir","Z","orb")[0]
        orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","I:p","pdir","Z","orb")[0]

        soci_perp.append([0,0,0,orb1,orb3])
        soci_perp.append([0,0,0,orb2,orb4])
        soci_par.append([0,0,0,orb2,orb5])
        soci_par.append([0,0,0,orb1,orb6])
        soci_par.append([0,0,0,orb3,orb6])
        soci_par.append([0,0,0,orb4,orb5])




#
# Construct "frames" of each hopping type here
#   Here "frame" means coordinate axes and not timestep.
# These are linearly independent vectors used to enforce rotation symmetry later
#

def xyzvec1(ddir):
  if ddir=="x":
    return np.array([1.0,0,0])
  elif ddir=="y":
    return np.array([0,1.0,0])
  elif ddir=="z":
    return np.array([0,0,1.0])

def xyzvec2(ddir):
  if ddir=="z":
    return np.array([1.0,0,0])
  elif ddir=="x":
    return np.array([0,1.0,0])
  elif ddir=="y":
    return np.array([0,0,1.0])


def pbp_frame(hopping): 
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec2(pdir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def pbs_frame(hopping): 
  return np.identity(3)

def ip_para_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec2(pdir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ip_perp_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  ddir1 = lut.get("orb",orb1,"dir")[0]
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec1(ddir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def pppi_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec1(pdir1)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ppsig_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec2(pdir)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def spsig_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  atom_typ1 = lut.get("atom",atom1,"atom_typ")[0]
  if atom_typ1=="Pb":
    pdir = lut.get("orb",orb2,"pdir")[0].lower()
  else:
    pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec2(pdir)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ppn_i_to_pb_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec1(pdir2)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)

  majoraxis = np.argmax(np.abs(axis2))
  axissign  = np.sign(axis2[majoraxis])
  axis1 = axis1*axissign

  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ppn_pb_to_i_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec1(pdir2)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)

  majoraxis = np.argmax(np.abs(axis2))
  axissign  = np.sign(axis2[majoraxis])
  axis1 = axis1*axissign

  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def spn_pb_to_i_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  atom_typ1 = lut.get("atom",atom1,"atom_typ")[0]
  if atom_typ1=="Pb":
    pdir = lut.get("orb",orb2,"pdir")[0].lower()
  else:
    pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def nnpbpb_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec2(pdir)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def nnii_cage_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  
  for neighbor1 in lut.get("atom",atom1,"ipb_shell_1"):
    atom3, ddx3, ddy3,ddz3 = neighbor1
    for neighbor2 in lut.get("atom",atom2,"ipb_shell_1"):
      atom4, ddx4, ddy4,ddz4 = neighbor2
      if atom3==atom4 and dx+ddx4==ddx3 and dy+ddy4==ddy3 and dz+ddz4==ddz3:
        atom_ref = atom3
        dx_ref = ddx3
        dy_ref = ddy3
        dz_ref = ddz3

  axis1 = gfc.posvec0(lut,atom_ref)-gfc.posvec0(lut,atom1)+np.array([NX*dx_ref,NY*dy_ref,NZ*dz_ref])
  axis1 = axis1/np.linalg.norm(axis1)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  axis3 = axis3/np.linalg.norm(axis3)
  axis1 = np.cross(axis2,axis3)
  axis1 = axis1/np.linalg.norm(axis1)
  return np.stack((axis1,axis2,axis3))

def nnii_sig_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec2(pdir)
  axis2 = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))


def socpb_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec1(pdir2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def soci_par_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  ddir1 = lut.get("orb",orb1,"dir")[0]
  ddir2 = lut.get("orb",orb2,"dir")[0]
  if ddir1==pdir1:
    axis1 = xyzvec1(pdir1)
    axis2 = xyzvec1(pdir2)
  elif ddir2==pdir2:
    axis1 = xyzvec1(pdir2)
    axis2 = xyzvec1(pdir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def soci_perp_frame(hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec1(pdir2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))




#reference frames
pbp_frame0 = pbp_frame(pbp[0])
pbs_frame0 = pbs_frame(pbs[0])
ip_para_frame0 = ip_para_frame(ip_para[0])
ip_perp_frame0 = ip_perp_frame(ip_perp[0])
pppi_frame0 = pppi_frame(pppi[0])
ppsig_frame0 = ppsig_frame(ppsig[0])
spsig_frame0 = spsig_frame(spsig[0])
ppn_i_to_pb_frame0 = ppn_i_to_pb_frame(ppn_i_to_pb[0])
ppn_pb_to_i_frame0 = ppn_pb_to_i_frame(ppn_pb_to_i[0])
spn_pb_to_i_frame0 = spn_pb_to_i_frame(spn_pb_to_i[0])
nnpbpb_frame0 = nnpbpb_frame(nnpbpb[0])
nnii_cage_frame0 = nnii_cage_frame(nnii_cage[0])
nnii_sig_frame0 = nnii_sig_frame(nnii_sig[0])
socpb_frame0 = socpb_frame(socpb[0]) 
soci_par_frame0 = soci_par_frame(soci_par[0]) 
soci_perp_frame0 = soci_perp_frame(soci_perp[0]) 

#
# define feature functions here
# featurefn(n, (dx,dy,dz,orb1,orb2)) = list of features for the n-th frame 
#

def geom_featurefn(n,hopping):
  #returns displacements to all atoms as features
  dx,dy,dz,orb1,orb2 = hopping
  atom = lut.get("orb",orb1,"atom")
  return np.reshape(coord_ij_sorted[n,atom,:,:], -1)


sortaxis =  np.array([2*np.pi,1/np.sqrt(2.0),np.e]) #this can be anything incommensurate with lattice and constant
def sortshell(atom1, shell,rotation):
  #sorts the atoms in a shell according to vector projection along sortaxis
  axis = np.matmul(np.linalg.inv(rotation),sortaxis)
  values = []
  r1 = gfc.posvec0(lut,atom1)
  for neighbor in shell:
    atom2,dx,dy,dz = neighbor
    r2 = gfc.posvec0(lut,atom2)+np.array([NX*dx,NY*dy,NZ*dz])
    values.append( np.dot(r2-r1,axis))
  values = np.array(values)
  shell = np.array(shell)
  return shell[np.argsort(values)]


def singlesite_featurefn(n,atom1, rotation):
  #returns bondvectors/distances in shells surrounding an atom
  #  atoms in the shell are rotated to the reference frame and sorted accoring to vector projection along sortaxis
  #  cartesian coordinates of each displacement vectors are rotated to the reference frame
  atom1_typ = lut.get("atom",atom1,"atom_typ")[0]
  ans = []


  ###
  global pb_pb_min
  global pb_i_min
  global i_i_min
  global i_pb_min
  global pb_pb_max
  global i_i_max
  global pb_i_max
  global i_pb_max
  rij_pbpb=[]
  rij_pbi=[]
  rij_ii=[]
  rij_ipb=[]


  if atom1_typ=="Pb":
    for ishell in range(mshells_pbpb): 
      neighbors = np.array(lut.get("atom", atom1, "pbpb_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        atom2,ddx,ddy,ddz = neighbor
        #print (atom2)
        #print (ddx,ddy,ddz)
        sij = coord[n,atom2,0,:]-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
        sij = np.matmul(rotation,sij)
        #print (sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        #if srij  < epsilon:
        #  continue
          #        sij = np.concatenate((sij_sorted*srij/(rij+epsilon),np.array([srij])))
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))
        if rij > pb_pb_max:
          pb_pb_max=rij
        if rij < pb_pb_min:
          pb_pb_min=rij
        rij_pbpb.append (rij)          


    for ishell in range(mshells_pbi): 
      neighbors = np.array(lut.get("atom", atom1, "pbi_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        atom2,ddx,ddy,ddz = neighbor
        sij = coord[n,atom2,0,:]-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        #if srij  < epsilon:
        #  continue
#        sij = np.concatenate((sij*srij/(rij+epsilon),np.array([srij])))
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))
        if rij > pb_i_max:
          pb_i_max=rij
        if rij < pb_i_min:
          pb_i_min=rij
        rij_pbi.append (rij)

    for ishell in range(mshells_pbma): 
      neighbors = np.array(lut.get("atom", atom1, "pbma_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        natom,ddx,ddy,ddz = neighbor
        ix = lut.get("atom",natom,"ix")[0]
        iy = lut.get("atom",natom,"iy")[0]
        iz = lut.get("atom",natom,"iz")[0]
        catom = lut.get("atom_typ","C","ix",ix,"iy",iy,"iz",iz,"atom")[0]
        sij = 0.5*(coord[n,natom,0,:]+coord[n,catom,0,:])-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))

    for ishell in range(mshells_pbma): 
      neighbors = np.array(lut.get("atom", atom1, "pbma_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        natom,ddx,ddy,ddz = neighbor
        ix = lut.get("atom",natom,"ix")[0]
        iy = lut.get("atom",natom,"iy")[0]
        iz = lut.get("atom",natom,"iz")[0]
        catom = lut.get("atom_typ","C","ix",ix,"iy",iy,"iz",iz,"atom")[0]
        sij = coord[n,natom,0,:]-coord[n,catom,0,:]
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))


  elif atom1_typ=="I1" or atom1_typ=="I2" or atom1_typ=="I3":
    for ishell in range(mshells_ipb): 
      neighbors = np.array(lut.get("atom", atom1, "ipb_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        atom2,ddx,ddy,ddz = neighbor
        sij = coord[n,atom2,0,:]-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        #if srij  < epsilon:
        #  continue
#        sij = np.concatenate((sij*srij/(rij+epsilon),np.array([srij])))
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))

        if rij > i_pb_max:
          i_pb_max=rij
        if rij < i_pb_min:
          i_pb_min=rij
        rij_ipb.append(rij)



    for ishell in range(mshells_ii): 
      neighbors = np.array(lut.get("atom", atom1, "ii_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        atom2,ddx,ddy,ddz = neighbor
        #print(atom1+1,atom2+1)
        sij = coord[n,atom2,0,:]-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        #if srij  < epsilon:
        #  continue
          #        sij = np.concatenate((sij*srij/(rij+epsilon),np.array([srij])))
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))
        if rij > i_i_max:
          i_i_max=rij
        if rij < i_i_min:
          i_i_min=rij
        rij_ii.append(rij)

    for ishell in range(mshells_ima): 
      neighbors = np.array(lut.get("atom", atom1, "ima_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        natom,ddx,ddy,ddz = neighbor
        ix = lut.get("atom",natom,"ix")[0]
        iy = lut.get("atom",natom,"iy")[0]
        iz = lut.get("atom",natom,"iz")[0]
        catom = lut.get("atom_typ","C","ix",ix,"iy",iy,"iz",iz,"atom")[0]
        sij = 0.5*(coord[n,natom,0,:]+coord[n,catom,0,:])-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))

    for ishell in range(mshells_ima): 
      neighbors = np.array(lut.get("atom", atom1, "ima_shell_%s" % (ishell+1)))
      neighbors_sorted = sortshell(atom1,neighbors,rotation)
      for neighbor in neighbors_sorted:
        natom,ddx,ddy,ddz = neighbor
        ix = lut.get("atom",natom,"ix")[0]
        iy = lut.get("atom",natom,"iy")[0]
        iz = lut.get("atom",natom,"iz")[0]
        catom = lut.get("atom_typ","C","ix",ix,"iy",iy,"iz",iz,"atom")[0]
        sij = coord[n,natom,0,:]-coord[n,catom,0,:]
        sij = np.matmul(rotation,sij)
        rij = np.linalg.norm(sij)
        srij=smooth_rij (rij)
        sij = np.concatenate((sij,np.array([srij])))
        ans = np.concatenate((ans,sij))

    rij_ipb=np.array(rij_ipb)
    rij_ii=np.array(rij_ii)    
    rij_ipb=rij_ipb.round(decimals=3)
    rij_ii=rij_ii.round(decimals=3)
    #print ('I-Pb',rij_ipb)
    #print ('I-I',rij_ii)        

  return ans

def shell_intersect_featurefn(n,atom1,atom2,dx,dy,dz,rotation):
  atom1_typ = lut.get("atom",atom1,"atom_typ")[0]
  atom2_typ = lut.get("atom",atom2,"atom_typ")[0]
  ans = []

  if atom1_typ=="Pb":
    shellnames1 = ["pbpb_shell_","pbi_shell_","pbma_shell_"]
    mshells_list1 = [mshells_pbpb,mshells_pbi,mshells_pbma]
  elif atom1_typ=="I1" or atom1_typ=="I2" or atom1_typ=="I3":
    shellnames1 = ["ipb_shell_","ii_shell_","ima_shell_"]
    mshells_list1 = [mshells_ipb,mshells_ii,mshells_ima]

  if atom2_typ=="Pb":
    shellnames2 = ["pbpb_shell_","pbi_shell_","pbma_shell_"]
    mshells_list2 = [mshells_pbpb,mshells_pbi,mshells_pbma]
  elif atom2_typ=="I1" or atom2_typ=="I2" or atom2_typ=="I3":
    shellnames2 = ["ipb_shell_","ii_shell_","ima_shell_"]
    mshells_list2 = [mshells_ipb,mshells_ii,mshells_ima]

  neighbors1 = []
  for kk in range(len(shellnames1)):
    shellname1 = shellnames1[kk]
    mshells_val = mshells_list1[kk]
    for ishell in range(mshells_val): 
      for neighbor in np.array(lut.get("atom", atom1, shellname1+str(ishell+1) )):
        atom3,ddx,ddy,ddz = neighbor
        if not (atom3==atom2 and ddx==dx and ddy==dy and ddz==dz):
          neighbors1.append(neighbor)
 
  neighbors2 = []
  for kk in range(len(shellnames2)):
    shellname2 = shellnames2[kk]
    mshells_val = mshells_list2[kk]
    for ishell in range(mshells_val): 
      for neighbor in np.array(lut.get("atom", atom2, shellname2+str(ishell+1) )):
        atom3,ddx,ddy,ddz = neighbor
        if not(atom3==atom1 and ddx+dx==0 and ddy+dy==0 and ddz+dz==0):
          neighbors2.append(neighbor)


  neighbors1_intersect = []
  neighbors2_intersect = []

  for neighbor1 in neighbors1:
    atom3,ddx3,ddy3,ddz3 = neighbor1
    for neighbor2 in neighbors2:
      atom4,ddx4,ddy4,ddz4 = neighbor2
      if atom3==atom4 and dx+ddx4==ddx3 and dy+ddy4==ddy3 and dz+ddz4==ddz3 :
        neighbors1_intersect.append(neighbor1)
        neighbors2_intersect.append(neighbor2)
          
  neighbors1_sorted = sortshell(atom1,neighbors1_intersect,rotation)
  neighbors2_sorted = sortshell(atom2,neighbors2_intersect,rotation)

  #this is the bond between atom1 and atom2 - always added to features regardless of mshells
  sij = coord[n,atom2,0,:]-coord[n,atom1,0,:] + np.array([dx*LX,dy*LY,dz*LZ])
  sij = np.matmul(rotation,sij)
  rij = np.linalg.norm(sij)
  srij=smooth_rij (rij)
  sij = np.concatenate((sij,np.array([srij])))
  ans = np.concatenate((ans,sij))

  for neighbor in neighbors1_sorted:
    atom3,ddx,ddy,ddz = neighbor
    atom3_typ = lut.get("atom",atom3,"atom_typ")[0]
    if atom3_typ!="N":
      sij = coord[n,atom3,0,:]-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
      sij = np.matmul(rotation,sij)
      rij = np.linalg.norm(sij)
      srij=smooth_rij (rij)
      sij = np.concatenate((sij,np.array([srij])))
      ans = np.concatenate((ans,sij))
    elif atom3_typ=="N":
      ix = lut.get("atom",atom3,"ix")[0]
      iy = lut.get("atom",atom3,"iy")[0]
      iz = lut.get("atom",atom3,"iz")[0]
      catom = lut.get("atom_typ","C","ix",ix,"iy",iy,"iz",iz,"atom")[0]
      #vector to C-N midpoint
      sij = 0.5*(coord[n,atom3,0,:]+coord[n,catom,0,:])-coord[n,atom1,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
      sij = np.matmul(rotation,sij)
      rij = np.linalg.norm(sij)
      srij=smooth_rij (rij)
      sij = np.concatenate((sij,np.array([srij])))
      ans = np.concatenate((ans,sij))
      #C-N vector
      sij = coord[n,atom3,0,:]-coord[n,catom,0,:]
      sij = np.matmul(rotation,sij)
      rij = np.linalg.norm(sij)
      srij=smooth_rij (rij)
      sij = np.concatenate((sij,np.array([srij])))
      ans = np.concatenate((ans,sij))


  for neighbor in neighbors2_sorted:
    atom3,ddx,ddy,ddz = neighbor
    atom3_typ = lut.get("atom",atom3,"atom_typ")[0]
    if atom3_typ!="N":
      sij = coord[n,atom3,0,:]-coord[n,atom2,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
      sij = np.matmul(rotation,sij)
      rij = np.linalg.norm(sij)
      srij=smooth_rij (rij)
      sij = np.concatenate((sij,np.array([srij])))
      ans = np.concatenate((ans,sij))
    elif atom3_typ=="N":
      ix = lut.get("atom",atom3,"ix")[0]
      iy = lut.get("atom",atom3,"iy")[0]
      iz = lut.get("atom",atom3,"iz")[0]
      catom = lut.get("atom_typ","C","ix",ix,"iy",iy,"iz",iz,"atom")[0]
      #vector to C-N midpoint
      sij = 0.5*(coord[n,atom3,0,:]+coord[n,catom,0,:])-coord[n,atom2,0,:] + np.array([ddx*LX,ddy*LY,ddz*LZ])
      sij = np.matmul(rotation,sij)
      rij = np.linalg.norm(sij)
      srij=smooth_rij (rij)
      sij = np.concatenate((sij,np.array([srij])))
      ans = np.concatenate((ans,sij))
      #C-N vector
      sij = coord[n,atom3,0,:]-coord[n,catom,0,:]
      sij = np.matmul(rotation,sij)
      rij = np.linalg.norm(sij)
      srij=smooth_rij (rij)
      sij = np.concatenate((sij,np.array([srij])))
      ans = np.concatenate((ans,sij))


  return ans



#geom feature functions below returns m-th nearest neighbor displacements as features
# Pb and I only

def pbp_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(pbp_frame0,np.linalg.inv(np.transpose(pbp_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)

def pbs_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(pbs_frame0,np.linalg.inv(np.transpose(pbs_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)

def ip_para_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(ip_para_frame0,np.linalg.inv(np.transpose(ip_para_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)

def ip_perp_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(ip_perp_frame0,np.linalg.inv(np.transpose(ip_perp_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)

def pppi_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(pppi_frame0,np.linalg.inv(np.transpose(pppi_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def ppsig_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(ppsig_frame0,np.linalg.inv(np.transpose(ppsig_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def spsig_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(spsig_frame0,np.linalg.inv(np.transpose(spsig_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def ppn_i_to_pb_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(ppn_i_to_pb_frame0,np.linalg.inv(np.transpose(ppn_i_to_pb_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def ppn_pb_to_i_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(ppn_pb_to_i_frame0,np.linalg.inv(np.transpose(ppn_pb_to_i_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def spn_pb_to_i_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(spn_pb_to_i_frame0,np.linalg.inv(np.transpose(spn_pb_to_i_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def nnpbpb_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(nnpbpb_frame0,np.linalg.inv(np.transpose(nnpbpb_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def nnii_cage_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(nnii_cage_frame0,np.linalg.inv(np.transpose(nnii_cage_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)

def nnii_sig_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom3 = lut.get("orb",orb3,"atom")[0]
  rotation = np.matmul(nnii_sig_frame0,np.linalg.inv(np.transpose(nnii_sig_frame(hopping))))
  #features1 = singlesite_featurefn(n,atom1,rotation)
  #features3 = singlesite_featurefn(n,atom3,rotation)
  #return np.concatenate((features1,features3))
  return shell_intersect_featurefn(n,atom1,atom3,dx,dy,dz,rotation)


def socpb_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  rotation = np.matmul(socpb_frame0,np.linalg.inv(np.transpose(socpb_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)

def soci_par_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  rotation = np.matmul(soci_par_frame0,np.linalg.inv(np.transpose(soci_par_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)

def soci_perp_geom_featurefn(n,hopping):
  dx,dy,dz,orb1,orb3 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  rotation = np.matmul(soci_perp_frame0,np.linalg.inv(np.transpose(soci_perp_frame(hopping))))
  return singlesite_featurefn(n,atom1,rotation)




#
# write feature files here
#

def write_features(hoppings,featurefn,filename):
  for n in range(nframe):
    ## becuase of the cuoff, number of features change with orbitlas. 
    # use minium value that works for all orbitals 
    nmin=10000000
    fout = open(filename+"_ts_"+str(n)+".dat",'w')
    for hopping in hoppings:
      dx = hopping[0]
      dy = hopping[1]
      dz = hopping[2]
      orb1 = hopping[3]
      orb2 = hopping[4]
      features = featurefn(n, hopping)
      if len(features) < nmin:
        nmin=len(features)
      fout.write("%s %s %s %s %s " % (dx, dy, dz, orb1, orb2))
      for feature in features:
        fout.write("%s " % feature)
      fout.write('\n')
    fout.close()
    #    print ('nmin', nmin)
    a=np.loadtxt(filename+"_ts_"+str(n)+".dat",usecols=[*range(0,nmin+5)])
    np.savetxt(filename+"_ts_"+str(n)+".dat",a)



write_features(pbp,pbp_geom_featurefn,"pbp")
write_features(pbs,pbs_geom_featurefn,"pbs")
write_features(ip_para,ip_para_geom_featurefn,"ip_para")
write_features(ip_perp,ip_perp_geom_featurefn,"ip_perp")

write_features(pppi,pppi_geom_featurefn,"pppi")
write_features(ppsig,ppsig_geom_featurefn,"ppsig")
write_features(spsig,spsig_geom_featurefn,"spsig")

write_features(ppn_i_to_pb,ppn_i_to_pb_geom_featurefn,"ppn_i_to_pb")
write_features(ppn_pb_to_i,ppn_pb_to_i_geom_featurefn,"ppn_pb_to_i")
write_features(spn_pb_to_i,spn_pb_to_i_geom_featurefn,"spn_pb_to_i")

#write_features(nnpbpb,nnpbpb_geom_featurefn,"nnpbpb")
#write_features(nnii_cage,nnii_cage_geom_featurefn,"nnii_cage")
#write_features(nnii_sig,nnii_sig_geom_featurefn,"nnii_sig")

write_features(socpb,socpb_geom_featurefn,"socpb")
write_features(soci_par,soci_par_geom_featurefn,"soci_par")
write_features(soci_perp,soci_perp_geom_featurefn,"soci_perp")


print (f" {pb_pb_min}  Pb-Pb  {pb_pb_max}")
print (f" {pb_i_min}  Pb-I  {pb_i_max}")
print (f" {i_i_min}  I-I  {i_i_max}")
print (f" {i_pb_min}  I-Pb  {i_pb_max}")
