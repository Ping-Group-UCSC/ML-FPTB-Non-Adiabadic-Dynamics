import numpy as np

# This script constructs shells of n-th neighbour atoms from a central Pb or I atom in an ideal perovskite structure

#search region:
mx = 7
my = 7
mz = 7

rpbpb = []
rpbi = []
rpbma = []

ri1pb = []
ri2pb = []
ri3pb = []

ri1i = []
ri2i = []
ri3i = []

ri1ma = []
ri2ma = []
ri3ma = []


rpbpb_idx = []
rpbi_idx = []
rpbma_idx = []

ri1pb_idx = []
ri2pb_idx = []
ri3pb_idx = []

ri1i_idx = []
ri2i_idx = []
ri3i_idx = []

ri1ma_idx = []
ri2ma_idx = []
ri3ma_idx = []

ri1_0 = np.array([0.5,0,0])
ri2_0 = np.array([0,0.5,0])
ri3_0 = np.array([0,0,0.5])


for ix in range(-mx,mx+1):
  for iy in range(-my,my+1):
    for iz in range(-mz,mz+1):
      rpb = np.array([ix,iy,iz])
      ri1 = np.array([ix+0.5,iy,iz])
      ri2 = np.array([ix,iy+0.5,iz])
      ri3 = np.array([ix,iy,iz+0.5])
      rma = np.array([ix+0.5,iy+0.5,iz+0.5])

      rpbpb.append(np.linalg.norm(rpb))

      rpbi.append(np.linalg.norm(ri1))
      rpbi.append(np.linalg.norm(ri2))
      rpbi.append(np.linalg.norm(ri3))

      rpbma.append(np.linalg.norm(rma))

      ri1pb.append(np.linalg.norm(rpb-ri1_0))
      ri2pb.append(np.linalg.norm(rpb-ri2_0))
      ri3pb.append(np.linalg.norm(rpb-ri3_0))

      ri1ma.append(np.linalg.norm(rma-ri1_0))
      ri2ma.append(np.linalg.norm(rma-ri2_0))
      ri3ma.append(np.linalg.norm(rma-ri3_0))

      ri1i.append(np.linalg.norm(ri1-ri1_0))
      ri1i.append(np.linalg.norm(ri2-ri1_0))
      ri1i.append(np.linalg.norm(ri3-ri1_0))

      ri2i.append(np.linalg.norm(ri1-ri2_0))
      ri2i.append(np.linalg.norm(ri2-ri2_0))
      ri2i.append(np.linalg.norm(ri3-ri2_0))

      ri3i.append(np.linalg.norm(ri1-ri3_0))
      ri3i.append(np.linalg.norm(ri2-ri3_0))
      ri3i.append(np.linalg.norm(ri3-ri3_0))


      rpbpb_idx.append([ix,iy,iz,0])

      rpbi_idx.append([ix,iy,iz,1])
      rpbi_idx.append([ix,iy,iz,2])
      rpbi_idx.append([ix,iy,iz,3])

      rpbma_idx.append([ix,iy,iz,0])

      ri1pb_idx.append([ix,iy,iz,0])
      ri2pb_idx.append([ix,iy,iz,0])
      ri3pb_idx.append([ix,iy,iz,0])

      ri1ma_idx.append([ix,iy,iz,0])
      ri2ma_idx.append([ix,iy,iz,0])
      ri3ma_idx.append([ix,iy,iz,0])

      ri1i_idx.append([ix,iy,iz,1])
      ri1i_idx.append([ix,iy,iz,2])
      ri1i_idx.append([ix,iy,iz,3])

      ri2i_idx.append([ix,iy,iz,1])
      ri2i_idx.append([ix,iy,iz,2])
      ri2i_idx.append([ix,iy,iz,3])

      ri3i_idx.append([ix,iy,iz,1])
      ri3i_idx.append([ix,iy,iz,2])
      ri3i_idx.append([ix,iy,iz,3])

rpbpb = np.array(rpbpb)
rpbi = np.array(rpbi)
rpbma = np.array(rpbma)

ri1pb = np.array(ri1pb)
ri2pb = np.array(ri2pb)
ri3pb = np.array(ri3pb)

ri1ma = np.array(ri1ma)
ri2ma = np.array(ri2ma)
ri3ma = np.array(ri3ma)

ri1i = np.array(ri1i)
ri2i = np.array(ri2i)
ri3i = np.array(ri3i)



rpbpb_idx = np.array(rpbpb_idx,dtype=int)
rpbi_idx = np.array(rpbi_idx,dtype=int)
rpbma_idx = np.array(rpbma_idx,dtype=int)

ri1pb_idx = np.array(ri1pb_idx,dtype=int)
ri2pb_idx = np.array(ri2pb_idx,dtype=int)
ri3pb_idx = np.array(ri3pb_idx,dtype=int)

ri1ma_idx = np.array(ri1ma_idx,dtype=int)
ri2ma_idx = np.array(ri2ma_idx,dtype=int)
ri3ma_idx = np.array(ri3ma_idx,dtype=int)

ri1i_idx = np.array(ri1i_idx,dtype=int)
ri2i_idx = np.array(ri2i_idx,dtype=int)
ri3i_idx = np.array(ri3i_idx,dtype=int)


def nshells(rs):
  eps = 1e-9
  ans = [0]

  for i in range(1,len(rs)):
    if rs[i] > rs[i-1]+eps:
      ans.append(i)
  return ans



np.savetxt("rpbpb_shells.dat",rpbpb[np.argsort(rpbpb)])
np.savetxt("rpbi_shells.dat",rpbi[np.argsort(rpbi)])
np.savetxt("rpbma_shells.dat",rpbma[np.argsort(rpbma)])

np.savetxt("ri1pb_shells.dat",ri1pb[np.argsort(ri1pb)])
np.savetxt("ri2pb_shells.dat",ri2pb[np.argsort(ri2pb)])
np.savetxt("ri3pb_shells.dat",ri3pb[np.argsort(ri3pb)])

np.savetxt("ri1ma_shells.dat",ri1ma[np.argsort(ri1ma)])
np.savetxt("ri2ma_shells.dat",ri2ma[np.argsort(ri2ma)])
np.savetxt("ri3ma_shells.dat",ri3ma[np.argsort(ri3ma)])

np.savetxt("ri1i_shells.dat",ri1i[np.argsort(ri1i)])
np.savetxt("ri2i_shells.dat",ri2i[np.argsort(ri2i)])
np.savetxt("ri3i_shells.dat",ri3i[np.argsort(ri3i)])


np.savetxt("n_pbpb_shells.dat",nshells(rpbpb[np.argsort(rpbpb)]),fmt="%8i")
np.savetxt("n_pbi_shells.dat",nshells(rpbi[np.argsort(rpbi)]),fmt="%8i")
np.savetxt("n_pbma_shells.dat",nshells(rpbma[np.argsort(rpbma)]),fmt="%8i")

np.savetxt("n_i1pb_shells.dat",nshells(ri1pb[np.argsort(ri1pb)]),fmt="%8i")
np.savetxt("n_i2pb_shells.dat",nshells(ri2pb[np.argsort(ri2pb)]),fmt="%8i")
np.savetxt("n_i3pb_shells.dat",nshells(ri3pb[np.argsort(ri3pb)]),fmt="%8i")

np.savetxt("n_i1ma_shells.dat",nshells(ri1ma[np.argsort(ri1ma)]),fmt="%8i")
np.savetxt("n_i2ma_shells.dat",nshells(ri2ma[np.argsort(ri2ma)]),fmt="%8i")
np.savetxt("n_i3ma_shells.dat",nshells(ri3ma[np.argsort(ri3ma)]),fmt="%8i")

np.savetxt("n_i1i_shells.dat",nshells(ri1i[np.argsort(ri1i)]),fmt="%8i")
np.savetxt("n_i2i_shells.dat",nshells(ri2i[np.argsort(ri2i)]),fmt="%8i")
np.savetxt("n_i3i_shells.dat",nshells(ri3i[np.argsort(ri3i)]),fmt="%8i")



np.savetxt("pbpb_shells.dat",rpbpb_idx[np.argsort(rpbpb),:],fmt="%4i")
np.savetxt("pbi_shells.dat",rpbi_idx[np.argsort(rpbi),:],fmt="%4i")
np.savetxt("pbma_shells.dat",rpbma_idx[np.argsort(rpbma),:],fmt="%4i")

np.savetxt("i1pb_shells.dat",ri1pb_idx[np.argsort(ri1pb),:],fmt="%4i")
np.savetxt("i2pb_shells.dat",ri2pb_idx[np.argsort(ri2pb),:],fmt="%4i")
np.savetxt("i3pb_shells.dat",ri3pb_idx[np.argsort(ri3pb),:],fmt="%4i")

np.savetxt("i1ma_shells.dat",ri1ma_idx[np.argsort(ri1ma),:],fmt="%4i")
np.savetxt("i2ma_shells.dat",ri2ma_idx[np.argsort(ri2ma),:],fmt="%4i")
np.savetxt("i3ma_shells.dat",ri3ma_idx[np.argsort(ri3ma),:],fmt="%4i")

np.savetxt("i1i_shells.dat",ri1i_idx[np.argsort(ri1i),:],fmt="%4i")
np.savetxt("i2i_shells.dat",ri2i_idx[np.argsort(ri2i),:],fmt="%4i")
np.savetxt("i3i_shells.dat",ri3i_idx[np.argsort(ri3i),:],fmt="%4i")


