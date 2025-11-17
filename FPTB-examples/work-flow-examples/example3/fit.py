import numpy as np
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"../../common")
import gen_features_common as gfc

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



#
# initialize lookuptables
#


lut = init.makelut("atom_index.dat","orb_index.dat")

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
# Initial guesses for parameters
#

vcoul_params_on_pb = (-0.2, -0.5, 0.08, 0.01, -11500., 45000., 47600., 171., -314., -517., 2691., -1950., -1710., -4.0)
vcoul_params_on_br = (-0.2, -0.04, 0.04, 0.04, -300., 2000., 1770., 410., 3500., 3600., 534.,-200., -300., -2.3)
vcoul_params_on_brp = (-0.2, -0.04, 0.04, 0.04, -300., 2000., 1770., 410., 3500., 3600., 534., -200., -300., -0.3, -2.3)
vcoul_params_soc = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)


lut.set("hoptype", "pbs", "params0", (-0.3, -0.004, -9500., 48600., -320., -3200, 7000., -4652., -10.))
lut.set("hoptype", "pbp", "params0", vcoul_params_on_pb)
lut.set("hoptype", "brp_par", "params0", vcoul_params_on_brp)
lut.set("hoptype", "brp_perp", "params0", vcoul_params_on_br)

lut.set("hoptype", "spsig", "params0", (1.0,-1.0,0.0,0.0))
lut.set("hoptype", "ppsig", "params0", (1.0,-1.0,0.0,0.0))
lut.set("hoptype", "pppi", "params0", (1.0,-1.0,0.0,0.0,0.0))

lut.set("hoptype", "ppn_br_to_pb", "params0", (1.0,-1.0,0.0))
lut.set("hoptype", "ppn_pb_to_br", "params0", (1.0,-1.0,0.0))
lut.set("hoptype", "spn_pb_to_br", "params0", (1.0,-1.0,0.0))

lut.set("hoptype", "socpb", "params0", vcoul_params_soc)
lut.set("hoptype", "socbr_perp", "params0", vcoul_params_soc)
lut.set("hoptype", "socbr_par", "params0", vcoul_params_soc)


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
# Collect data
#

timesteps = ["100fsmd", "200fsmd"]

#pbs_data = []
#pbp_data = []
#brp_par_data = []
#brp_par_data = []
#spsig_data = []
#ppsig_data = []
#pppi_data = []
#socpb_data = []
#socbr_par_data = []
#socbr_perp_data = []
#
#pbs_features = []
#pbp_features = []
#brp_par_features = []
#brp_par_features = []
#spsig_features = []
#ppsig_features = []
#pppi_features = []
#socpb_features = []
#socbr_par_features = []
#socbr_perp_features = []

for timestep in timesteps:
  print("reading " + timestep)
  coord = np.loadtxt("wanndata/"+timestep+"/atompos.dat")
  coord = coord * np.array([LX, LY, LZ])
  #translate atoms to be centered on standard positions:
  for ix in range(NX):
    for iy in range(NY):
      for iz in range(NZ):
        for atom in lut.get("ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
          coord[atom,:] = stdpos.translate_posvec0(lut, atom, coord[atom,:], LX0, LY0, LZ0)

  #print(coord)


  hrdat = open("wanndata/"+timestep+"/test_hr.dat","r")
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

    hopping = (dx,dy,dz,orb1,orb2)

    match = lut.get("hopping",hopping,"hoptype")
    if len(match)==0:
      continue
    hoptype = match[0]

    
    phase_match = lut.get("hopping",hopping,"phase")
    if len(phase_match)==0:
      continue
    hval = np.real(lut.get("hopping",hopping,"phase")[0] * (reh+1.0j*imh) )

    #print(hopping)

    lut.set("hoptype",hoptype,"hdata",hval)
    featurefn = lut.get("hoptype",hoptype,"featurefn")[0]
    feats = featurefn(lut,NX,NY,NZ,cell,qs,hopping,coord)
    lut.set("hoptype",hoptype,"featdata",feats)



  hrdat.close()

#
#do fit
#

hoptypes = ["pbs","pbp","brp_par","brp_perp","spsig","ppsig","pppi","ppn_br_to_pb", "ppn_pb_to_br" ,"spn_pb_to_br", "socpb", "socbr_par","socbr_perp"]
#hoptypes = []

for hoptype in hoptypes:
  print("fitting " + hoptype)
  hdata = lut.get("hoptype",hoptype,"hdata")
  featdata = lut.get("hoptype",hoptype,"featdata")
  featdata = np.transpose( np.array(featdata))
  print(featdata)
  model = lut.get("hoptype",hoptype,"model")[0]
  p0 = lut.get("hoptype",hoptype,"params0")[0]
  fitparams, fitcov = curve_fit(model,featdata,hdata,p0)

  pred_hdata = model(featdata, *fitparams)


  np.savetxt(hoptype+"_fit.dat",np.stack((hdata,pred_hdata),axis=1))
  np.savetxt(hoptype+"_params.dat",fitparams)

  plt.clf()
  plt.plot(hdata,pred_hdata,'o')
  plt.title(hoptype)
  plt.xlabel("Wannier (eV)")
  plt.ylabel("Model (eV)")
  plt.savefig(hoptype+"_fit.png")
