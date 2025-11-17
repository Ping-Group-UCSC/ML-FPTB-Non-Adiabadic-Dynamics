import numpy as np
from vsr_derivatives import *
import sys
sys.path.insert(1,"../../common")
import gen_features_common as gfc
import standard_positions as stdpos

#initialize charges

def make_charges(lut,NN):
  qs = []
  for atom in range(NN):
    atom_typ = lut.get("atom",atom, "atom_typ")[0]
    if atom_typ == "Pb":
      qs.append(2.0)
    elif atom_typ=="Br1" or atom_typ=="Br2" or atom_typ=="Br3" :
      qs.append(-1.0)
    elif atom_typ=="Cs":
      qs.append(1.0)
  return qs



#
# Construct "frames" of each hopping type here
#   Here "frame" means coordinate axes 
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


def pbp_frame(lut,NX,NY,NZ,hopping): 
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec2(pdir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def pbs_frame(lut,NX,NY,NZ,hopping): 
  return np.identity(3)

def brp_par_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec2(pdir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def brp_perp_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  ddir1 = lut.get("orb",orb1,"dir")[0]
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec1(ddir1)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def pppi_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec1(pdir1)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ppsig_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec2(pdir)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def spsig_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  atom_typ1 = lut.get("atom",atom1,"atom_typ")[0]
  if atom_typ1=="Pb":
    pdir = lut.get("orb",orb2,"pdir")[0].lower()
  else:
    pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec2(pdir)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ppn_br_to_pb_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec1(pdir2)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)

  majoraxis = np.argmax(np.abs(axis2))
  axissign  = np.sign(axis2[majoraxis])
  axis1 = axis1*axissign

  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def ppn_pb_to_br_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  axis1 = xyzvec1(pdir2)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)

  majoraxis = np.argmax(np.abs(axis2))
  axissign  = np.sign(axis2[majoraxis])
  axis1 = axis1*axissign

  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def spn_pb_to_br_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  atom_typ1 = lut.get("atom",atom1,"atom_typ")[0]
  if atom_typ1=="Pb":
    pdir = lut.get("orb",orb2,"pdir")[0].lower()
  else:
    pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec1(pdir)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def nnpbpb_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec2(pdir)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def nnii_cage_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  
  for neighbor1 in lut.get("atom",atom1,"brpb_shell_1"):
    atom3, ddx3, ddy3,ddz3 = neighbor1
    for neighbor2 in lut.get("atom",atom2,"brpb_shell_1"):
      atom4, ddx4, ddy4,ddz4 = neighbor2
      if atom3==atom4 and dx+ddx4==ddx3 and dy+ddy4==ddy3 and dz+ddz4==ddz3:
        atom_ref = atom3
        dx_ref = ddx3
        dy_ref = ddy3
        dz_ref = ddz3

  axis1 = stdpos.posvec0(lut,atom_ref)-stdpos.posvec0(lut,atom1)+np.array([NX*dx_ref,NY*dy_ref,NZ*dz_ref])
  axis1 = axis1/np.linalg.norm(axis1)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  axis3 = axis3/np.linalg.norm(axis3)
  axis1 = np.cross(axis2,axis3)
  axis1 = axis1/np.linalg.norm(axis1)
  return np.stack((axis1,axis2,axis3))

def nnii_sig_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  pdir = lut.get("orb",orb1,"pdir")[0].lower()
  axis1 = xyzvec2(pdir)
  axis2 = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
  axis2 = axis2/np.linalg.norm(axis2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))


def socpb_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec1(pdir2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def socbr_par_frame(lut,NX,NY,NZ,hopping):
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

def socbr_perp_frame(lut,NX,NY,NZ,hopping):
  dx,dy,dz,orb1,orb2 = hopping
  pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
  pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  axis1 = xyzvec1(pdir1)
  axis2 = xyzvec1(pdir2)
  axis3 = np.cross(axis1,axis2)
  return np.stack((axis1,axis2,axis3))

def set_frames(lut):
  lut.set("hoptype", "pbp", "frame",pbp_frame)
  lut.set("hoptype", "pbs", "frame",pbs_frame)
  lut.set("hoptype", "brp_par", "frame", brp_par_frame)
  lut.set("hoptype", "brp_perp", "frame", brp_perp_frame)
  lut.set("hoptype", "pppi", "frame", pppi_frame)
  lut.set("hoptype", "ppsig", "frame", ppsig_frame)
  lut.set("hoptype", "spsig", "frame", spsig_frame)
  lut.set("hoptype", "ppn_br_to_pb", "frame",ppn_br_to_pb_frame)
  lut.set("hoptype", "ppn_pb_to_br", "frame",ppn_pb_to_br_frame)
  lut.set("hoptype", "spn_pb_to_br", "frame",spn_pb_to_br_frame)
  lut.set("hoptype", "socpb", "frame",socpb_frame)
  lut.set("hoptype", "socbr_par", "frame",socbr_par_frame)
  lut.set("hoptype", "socbr_perp", "frame",socbr_perp_frame)


#reference frames
#pbp_frame0 = pbp_frame(pbp[0])
#pbs_frame0 = pbs_frame(pbs[0])
#ip_para_frame0 = ip_para_frame(ip_para[0])
#ip_perp_frame0 = ip_perp_frame(ip_perp[0])
#pppi_frame0 = pppi_frame(pppi[0])
#ppsig_frame0 = ppsig_frame(ppsig[0])
#spsig_frame0 = spsig_frame(spsig[0])
#ppn_i_to_pb_frame0 = ppn_i_to_pb_frame(ppn_i_to_pb[0])
#ppn_pb_to_i_frame0 = ppn_pb_to_i_frame(ppn_pb_to_i[0])
#spn_pb_to_i_frame0 = spn_pb_to_i_frame(spn_pb_to_i[0])
#nnpbpb_frame0 = nnpbpb_frame(nnpbpb[0])
#nnii_cage_frame0 = nnii_cage_frame(nnii_cage[0])
#nnii_sig_frame0 = nnii_sig_frame(nnii_sig[0])
#socpb_frame0 = socpb_frame(socpb[0]) 
#soci_par_frame0 = soci_par_frame(soci_par[0]) 
#soci_perp_frame0 = soci_perp_frame(soci_perp[0]) 



def vcoul_vsr(lut,cell,qs,hopping,atompos):
  #assumes that atompos has been brought into standard positons
  orb1 = hopping[3]
  atom1 = lut.get("orb",orb1,"atom")[0]
<<<<<<< HEAD
  vcoul = gfc.ewald_phi(atom1,cell,qs,np.transpose(atompos),12.0,2,2)
  return (vcoul,)
=======
  vcoul = gfc.ewald_phi(atom1,cell,qs,np.transpose(atompos),48,2,2)

  qs = np.array(qs)
  Brmsk = qs == -1.0
  Pbmsk = qs == 2.0
  Csmsk = qs == 1.0
  
  vsrPb = get_vsr(atom1,cell,qs, np.transpose(atompos), Pbmsk, 1.33, 2)
  vsrBr = get_vsr(atom1,cell,qs, np.transpose(atompos), Brmsk, 1.82, 2)
  vsrCs = get_vsr(atom1,cell,qs, np.transpose(atompos), Csmsk, 1.81, 2)

  return (vcoul, vsrPb, vsrBr, vsrCs)

def pbs_coul_vsr_featurefn(lut,NX,NY,NZ,cell,qs,hopping,atompos):

  orb1 = hopping[3]
  atom1 = lut.get("orb",orb1,"atom")[0]

  qs = np.array(qs)
  Brmsk = qs == -1.0
  Pbmsk = qs == 2.0
  Csmsk = qs == 1.0

  vcoul, vsrPb, vsrBr, vsrCs = vcoul_vsr(lut,cell,qs,hopping,atompos)
  vcould4x, vcould4y, vcould4z = get_d4coul(atom1, cell, qs, np.transpose(atompos), 2)
  vsrPbd2x, vsrPbd2y, vsrPbd2z = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Pbmsk, 1.33,2)
  vsrBrd2x, vsrBrd2y, vsrBrd2z = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Brmsk, 1.82,2)
  vsrCsd2x, vsrCsd2y, vsrCsd2z = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Csmsk, 1.81,2)

  vcould4 = vcould4x + vcould4y + vcould4z
  vsrPbd2 = vsrPbd2x + vsrPbd2y + vsrPbd2z
  vsrBrd2 = vsrBrd2x + vsrBrd2y + vsrBrd2z
  vsrCsd2 = vsrCsd2x + vsrCsd2y + vsrCsd2z

  return (vcoul, vcould4, vsrPb, vsrPbd2, vsrBr, vsrBrd2, vsrCs, vsrCsd2) 

def pbp_coul_vsr_featurefn(lut,NX,NY,NZ,cell,qs,hopping,atompos):

  orb1 = hopping[3]
  atom1 = lut.get("orb",orb1,"atom")[0]

  qs = np.array(qs)
  Brmsk = qs == -1.0
  Pbmsk = qs == 2.0
  Csmsk = qs == 1.0

  direction = pbp_frame(lut,NX,NY,NZ,hopping)[0]
  dmsk = np.array(direction, dtype=bool)
  vcoul, vsrPb, vsrBr, vsrCs = vcoul_vsr(lut,cell,qs,hopping,atompos)
  vcould2r = get_d2coul(atom1, cell, qs, np.transpose(atompos), 2)
  vcould4r = get_d4coul(atom1, cell, qs, np.transpose(atompos), 2)
  vsrPbd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Pbmsk, 1.33,2)
  vsrBrd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Brmsk, 1.82,2)
  vsrCsd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Csmsk, 1.81,2)

  vcould2x = vcould2r[dmsk][0]
  vcould4x = vcould4r[dmsk][0]
  vsrPbd2x = vsrPbd2r[dmsk][0]
  vsrBrd2x = vsrBrd2r[dmsk][0]
  vsrCsd2x = vsrCsd2r[dmsk][0]

  vcould4yz = np.sum(vcould4r[~dmsk])
  vsrPbd2yz = np.sum(vsrPbd2r[~dmsk])
  vsrBrd2yz = np.sum(vsrBrd2r[~dmsk])
  vsrCsd2yz = np.sum(vsrCsd2r[~dmsk])

  return (vcoul, vcould2x, vcould4x, vcould4yz, vsrPb, vsrPbd2x, vsrPbd2yz,\
         vsrBr, vsrBrd2x, vsrBrd2yz, vsrCs, vsrCsd2x, vsrCsd2yz)

def brp_par_coul_vsr_featurefn(lut,NX,NY,NZ,cell,qs,hopping,atompos):

  orb1 = hopping[3]
  atom1 = lut.get("orb",orb1,"atom")[0]

  qs = np.array(qs)
  Brmsk = qs == -1.0
  Pbmsk = qs == 2.0
  Csmsk = qs == 1.0

  direction = brp_par_frame(lut,NX,NY,NZ,hopping)[0]
  dmsk = np.array(direction, dtype=bool)

  vcoul, vsrPb, vsrBr, vsrCs = vcoul_vsr(lut,cell,qs,hopping,atompos)
  vcould2r = get_d2coul(atom1, cell, qs, np.transpose(atompos), 2)
  vcould4r = get_d4coul(atom1, cell, qs, np.transpose(atompos), 2)
  vsrPbd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Pbmsk, 1.33,2)
  vsrBrd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Brmsk, 1.82,2)
  vsrCsd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Csmsk, 1.81,2)

  vcould2x = vcould2r[dmsk][0]
  vcould4x = vcould4r[dmsk][0]
  vsrPbd2x = vsrPbd2r[dmsk][0]
  vsrBrd2x = vsrBrd2r[dmsk][0]
  vsrCsd2x = vsrCsd2r[dmsk][0]

  vcould4yz = np.sum(vcould4r[~dmsk])
  vsrPbd2yz = np.sum(vsrPbd2r[~dmsk])
  vsrBrd2yz = np.sum(vsrBrd2r[~dmsk])
  vsrCsd2yz = np.sum(vsrCsd2r[~dmsk])

  return vcoul, vcould2x, vcould4x, vcould4yz, vsrPb, vsrPbd2x, vsrPbd2yz,\
         vsrBr, vsrBrd2x, vsrBrd2yz, vsrCs, vsrCsd2x, vsrCsd2yz

def brp_perp_coul_vsr_featurefn(lut,NX,NY,NZ,cell,qs,hopping,atompos):

  orb1 = hopping[3]
  atom1 = lut.get("orb",orb1,"atom")[0]

  qs = np.array(qs)
  Brmsk = qs == -1.0
  Pbmsk = qs == 2.0
  Csmsk = qs == 1.0

  direction = brp_perp_frame(lut,NX,NY,NZ,hopping)[0]
  dmsk = np.array(direction, dtype=bool)

  vcoul, vsrPb, vsrBr, vsrCs = vcoul_vsr(lut,cell,qs,hopping,atompos)
  vcould2r = get_d2coul(atom1, cell, qs, np.transpose(atompos), 2)
  vcould4r = get_d4coul(atom1, cell, qs, np.transpose(atompos), 2)
  vsrPbd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Pbmsk, 1.33,2)
  vsrBrd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Brmsk, 1.82,2)
  vsrCsd2r = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Csmsk, 1.81,2)

  vcould2x = vcould2r[dmsk][0]
  vcould4x = vcould4r[dmsk][0]
  vsrPbd2x = vsrPbd2r[dmsk][0]
  vsrBrd2x = vsrBrd2r[dmsk][0]
  vsrCsd2x = vsrCsd2r[dmsk][0]

  vcould4yz = np.sum(vcould4r[~dmsk])
  vsrPbd2yz = np.sum(vsrPbd2r[~dmsk])
  vsrBrd2yz = np.sum(vsrBrd2r[~dmsk])
  vsrCsd2yz = np.sum(vsrCsd2r[~dmsk])

  return vcoul, vcould2x, vcould4x, vcould4yz, vsrPb, vsrPbd2x, vsrPbd2yz,\
         vsrBr, vsrBrd2x, vsrBrd2yz, vsrCs, vsrCsd2x, vsrCsd2yz

def soc_coul_vsr_featurefn(lut,NX,NY,NZ,cell,qs,hopping,atompos):

  orb1 = hopping[3]
  atom1 = lut.get("orb",orb1,"atom")[0]

  qs = np.array(qs)
  Brmsk = qs == -1.0
  Pbmsk = qs == 2.0
  Csmsk = qs == 1.0

  vcoul, vsrPb, vsrBr, vsrCs = vcoul_vsr(lut,cell,qs,hopping,atompos)
  vcould2x, vcould2y, vcould2z = get_d2coul(atom1, cell, qs, np.transpose(atompos), 2)
  vcould4x, vcould4y, vcould4z = get_d4coul(atom1, cell, qs, np.transpose(atompos), 2)
  vsrPbd2x, vsrPbd2y, vsrPbd2z = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Pbmsk, 1.33,2)
  vsrBrd2x, vsrBrd2y, vsrBrd2z = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Brmsk, 1.82,2)
  vsrCsd2x, vsrCsd2y, vsrCsd2z = get_d2vsr(atom1, cell, qs, np.transpose(atompos), Csmsk, 1.81,2)


  return vcoul, vcould2x, vcould2y, vcould4x, vcould4y, vcould4z, vsrPb, vsrPbd2x, vsrPbd2y,\
         vsrPbd2z, vsrBr, vsrBrd2x, vsrBrd2y, vsrBrd2z, vsrCs, vsrCsd2x, vsrCsd2y, vsrCsd2z
>>>>>>> 0fb0ef5a31035d6df026a2b6eea4e2e113dedf33

def xyz_featurefn(lut,NX,NY,NZ,cell,qs,hopping,atompos):
  #assumes that atompos has been brought into standard positons
  hoptype =lut.get("hopping", hopping,"hoptype")[0]
  framefn = lut.get("hoptype",hoptype,"frame")[0]
  fr = framefn(lut,NX,NY,NZ,hopping)
  dx,dy,dz,orb1,orb2 = hopping
  atom1 = lut.get("orb",orb1,"atom")[0]
  atom2 = lut.get("orb",orb2,"atom")[0]
  #if(hoptype=="spsig"):
  #  print(atompos[atom2,:])
  #  print(atompos[atom1,:])
  #  print(np.matmul(cell,np.array([dx,dy,dz])))
  rdiff = atompos[atom2,:] - atompos[atom1,:] + np.matmul(cell,np.array([dx,dy,dz]))
  xdiff = np.dot(rdiff,fr[1,:]) #axis 2 was defined in the frames to be the bonding direction
  ydiff = np.dot(rdiff,fr[0,:])
  zdiff = np.dot(rdiff,fr[2,:])

  return (xdiff,ydiff,zdiff)

def set_featurefns(lut):
  lut.set("hoptype", "pbs", "featurefn", pbs_coul_vsr_featurefn)
  lut.set("hoptype", "pbp", "featurefn", pbp_coul_vsr_featurefn)
  lut.set("hoptype", "brp_par", "featurefn", brp_par_coul_vsr_featurefn)
  lut.set("hoptype", "brp_perp", "featurefn", brp_perp_coul_vsr_featurefn)
  
  lut.set("hoptype", "spsig", "featurefn", xyz_featurefn)
  lut.set("hoptype", "ppsig", "featurefn", xyz_featurefn)
  lut.set("hoptype", "pppi", "featurefn", xyz_featurefn)
  
  lut.set("hoptype", "ppn_br_to_pb", "featurefn", xyz_featurefn)
  lut.set("hoptype", "ppn_pb_to_br", "featurefn", xyz_featurefn)
  lut.set("hoptype", "spn_pb_to_br", "featurefn", xyz_featurefn)
  
  lut.set("hoptype", "socpb", "featurefn", soc_coul_vsr_featurefn)
  lut.set("hoptype", "socbr_perp", "featurefn", soc_coul_vsr_featurefn)
  lut.set("hoptype", "socbr_par", "featurefn", soc_coul_vsr_featurefn)



