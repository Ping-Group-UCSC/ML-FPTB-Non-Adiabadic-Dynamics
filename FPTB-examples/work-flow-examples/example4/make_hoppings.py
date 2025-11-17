import numpy as np

import sys
sys.path.insert(1,"../../common")
import gen_features_common as gfc

sys.path.insert(1,"../example3")
import init
import load_shells
import hoppingtypes as ht
import models 
import featurefns
import phases
import standard_positions as stdpos

NN = 40 # 2x2x2 supercell, all atoms
NX = 2
NY = 2
NZ = 2
LX = 12.00445150246208 #supercell dimensions
LY = 12.00445150246208
LZ = 12.00445150246208
LX0 = LX/NX #unitcell dimensions
LY0 = LY/NY
LZ0 = LZ/NZ

cell = np.array([LX,LY,LZ]) * np.identity(3)

lut = init.makelut("../example3/atom_index.dat","../example3/orb_index.dat")

coord = np.loadtxt("../example3/wanndata/100fsmd/atompos.dat")
coord = coord * np.array([LX, LY, LZ])
#translate atoms to be centered on standard positions:
for ix in range(NX):
  for iy in range(NY):
    for iz in range(NZ):
      for atom in lut.get("ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
        coord[atom,:] = stdpos.translate_posvec0(lut, atom, coord[atom,:], LX0, LY0, LZ0)



#
#load shells of n-th neighbor atoms

load_shells.set_shells(lut,NX,NY,NZ)

# construct hoppings of each type here
# dx, dy, dz, orb1, orb2

ht.set_hoppings(lut,NX,NY,NZ)

#
# Construct "frames" of each hopping type here
#   Here "frame" means coordinate axes 
# These are linearly independent vectors used to enforce rotation symmetry later
#

featurefns.set_frames(lut)

#
# model functions
#

models.set_models(lut)

#
# Phase transformations to be applied to hopping parameters for fitting
#

phases.set_phases(lut,NX,NY,NZ)

#
# feature functions
#

qs = featurefns.make_charges(lut,NN)
featurefns.set_featurefns(lut)

#
# define a function to get the dn-spin hoppings
#
def dnspinhops(lut,hopping):

  dx = hopping[0]
  dy = hopping[1]
  dz = hopping[2]
  orb1 = hopping[3]
  orb2 = hopping[4]
  spin1 = lut.get("orb",orb1,"spin")[0]
  spin2 = lut.get("orb",orb2,"spin")[0]

  orb_typ1 = lut.get("orb",orb1,"orb_typ")[0]
  ix1 = lut.get("orb",orb1,"ix")[0]
  iy1 = lut.get("orb",orb1,"iy")[0]
  iz1 = lut.get("orb",orb1,"iz")[0]
  if orb_typ1=="Pb:s":
    orb3 = lut.get("orb_typ",orb_typ1,"ix",ix1,"iy",iy1,"iz",iz1,"spin","dn","orb")[0]
  if orb_typ1=="Pb:p":
    pdir1 = lut.get("orb",orb1,"pdir")[0]
    orb3 = lut.get("orb_typ",orb_typ1,"ix",ix1,"iy",iy1,"iz",iz1,"pdir",pdir1,"spin","dn","orb")[0]
  if orb_typ1=="Br:p":
    dir1 = lut.get("orb",orb1,"dir")[0]
    pdir1 = lut.get("orb",orb1,"pdir")[0]
    orb3 = lut.get("orb_typ",orb_typ1,"ix",ix1,"iy",iy1,"iz",iz1,"pdir",pdir1,"dir",dir1,"spin","dn","orb")[0]

  orb_typ2 = lut.get("orb",orb2,"orb_typ")[0]
  ix2 = lut.get("orb",orb2,"ix")[0]
  iy2 = lut.get("orb",orb2,"iy")[0]
  iz2 = lut.get("orb",orb2,"iz")[0]
  if orb_typ2=="Pb:s":
    orb4 = lut.get("orb_typ",orb_typ2,"ix",ix2,"iy",iy2,"iz",iz2,"spin","dn","orb")[0]
  if orb_typ2=="Pb:p":
    pdir2 = lut.get("orb",orb2,"pdir")[0]
    orb4 = lut.get("orb_typ",orb_typ2,"ix",ix2,"iy",iy2,"iz",iz2,"pdir",pdir2,"spin","dn","orb")[0]
  if orb_typ2=="Br:p":
    dir2 = lut.get("orb",orb2,"dir")[0]
    pdir2 = lut.get("orb",orb2,"pdir")[0]
    orb4 = lut.get("orb_typ",orb_typ2,"ix",ix2,"iy",iy2,"iz",iz2,"pdir",pdir2,"dir",dir2,"spin","dn","orb")[0]


  return [dx,dy,dz,orb3,orb4]



hoptypes = ["pbs","pbp","brp_par","brp_perp","spsig","ppsig","pppi","ppn_br_to_pb", "ppn_pb_to_br" ,"spn_pb_to_br", "socpb", "socbr_par","socbr_perp"]
all_hoppings = []

for hoptype in hoptypes:
  params = np.loadtxt("../example3/"+hoptype+"_params.dat")
  model = lut.get("hoptype",hoptype,"model")[0]

  featurefn = lut.get("hoptype",hoptype,"featurefn")[0]

  for hopping in lut.get("hoptype",hoptype,"hopping"):
    feats = featurefn(lut,NX,NY,NZ,cell,qs,hopping,coord)
    lut.set("hoptype",hoptype,"featdata",feats)


  featdata = lut.get("hoptype",hoptype,"featdata")
  featdata = np.transpose( np.array(featdata))

  pred_hdata = model(featdata, *params)

  i=0
  for hopping in lut.get("hoptype",hoptype,"hopping"):
    hval = pred_hdata[i] / lut.get("hopping",hopping,"phase")[0] 

    reh = np.real(hval )
    imh = np.imag(hval )
    i=i+1
    all_hoppings.append( np.concatenate((np.array(hopping),np.array([reh,imh]))))

    #put in the dn spins
    dx = hopping[0]
    dy = hopping[1]
    dz = hopping[2]
    orb1 = hopping[3]
    orb2 = hopping[4]
    spin1 = lut.get("orb",orb1,"spin")[0]
    spin2 = lut.get("orb",orb2,"spin")[0]

    if (spin1=="up" and spin2=="up" and hoptype!="socbr_par"  and hoptype!="socbr_perp" and  hoptype!="socpb"):
      all_hoppings.append( np.concatenate((np.array(dnspinhops(lut,hopping)),np.array([reh,imh]))))

np.savetxt("hopping.dat",all_hoppings)

#original hoppings

orig_hoppings = []
hrdat = open("../example3/wanndata/100fsmd/test_hr.dat","r")
for line in hrdat:
  ll = line.split()

  if len(ll) != 7:
    continue

  dx = int(ll[0])
  dy = int(ll[1])
  dz = int(ll[2])
  orb1 = int(ll[3])
  orb2 = int(ll[4])
  reh = float(ll[5])
  imh = float(ll[6])
  spin1 = lut.get("orb",orb1,"spin")[0]
  spin2 = lut.get("orb",orb2,"spin")[0]

  hopping = (dx,dy,dz,orb1,orb2)

  match = lut.get("hopping",hopping,"hoptype")
  if len(match)==0:
    continue
  hoptype = match[0]

  if hoptype in hoptypes:
    orig_hoppings.append([dx,dy,dz,orb1,orb2,reh,imh])

    if (spin1=="up" and spin2=="up" and hoptype!="socbr_par"  and hoptype!="socbr_perp" and  hoptype!="socpb"):
      orig_hoppings.append( np.concatenate((np.array(dnspinhops(lut,hopping)),np.array([reh,imh]))))

dnspinhops(lut,hopping)

np.savetxt("orig_hopping.dat",orig_hoppings)
