import numpy as np
#
# model functions
#

def pbs_model(features,*params):
  vcoul = features[0]
  vcould4 = features[1]

  vsrPb = features[2]
  vsrPbd2 = features[3]

  vsrBr = features[4]
  vsrBrd2 = features[5]

  vsrCs = features[6]
  vsrCsd2 = features[7]

  return params[0]*vcoul + params[1]*vcould4 + params[2]*vsrPb + params[3]*vsrPbd2 \
                         + params[4]*vsrBr + params[5]*vsrBrd2 \
                         + params[6]*vsrCs + params[7]*vsrCsd2 + params[8]


def pbp_model(features,*params):
  vcoul = features[0]
  vcould2x = features[1]
  vcould4x = features[2]
  vcould4yz = features[3]

  vsrPb = features[4]
  vsrPbd2x = features[5]
  vsrPbd2yz = features[6]

  vsrBr = features[7]
  vsrBrd2x = features[8]
  vsrBrd2yz = features[9]

  vsrCs = features[10]
  vsrCsd2x = features[11]
  vsrCsd2yz = features[12]

  return params[0]*vcoul + params[1]*vcould2x + params[2]*vcould4x + params[3]*vcould4yz \
                         + params[4]*vsrPb + params[5]*vsrPbd2x + params[6]*vsrPbd2yz \
                         + params[7]*vsrBr + params[8]*vsrBrd2x + params[9]*vsrBrd2yz \
                         + params[10]*vsrCs + params[11]*vsrCsd2x + params[12]*vsrCsd2yz + params[13]

def brp_perp_model(features,*params):
  vcoul = features[0]
  vcould2x = features[1]
  vcould4x = features[2]
  vcould4yz = features[3]

  vsrPb = features[4]
  vsrPbd2x = features[5]
  vsrPbd2yz = features[6]

  vsrBr = features[7]
  vsrBrd2x = features[8]
  vsrBrd2yz = features[9]

  vsrCs = features[10]
  vsrCsd2x = features[11]
  vsrCsd2yz = features[12]

  return params[0]*vcoul + params[1]*vcould2x + params[2]*vcould4x + params[3]*vcould4yz \
                         + params[4]*vsrPb + params[5]*vsrPbd2x + params[6]*vsrPbd2yz \
                         + params[7]*vsrBr + params[8]*vsrBrd2x + params[9]*vsrBrd2yz \
                         + params[10]*vsrCs + params[11]*vsrCsd2x + params[12]*vsrCsd2yz + params[13]

def brp_par_model(features,*params):
  vcoul = features[0]
  vcould2x = features[1]
  vcould4x = features[2]
  vcould4yz = features[3]

  vsrPb = features[4]
  vsrPbd2x = features[5]
  vsrPbd2yz = features[6]

  vsrBr = features[7]
  vsrBrd2x = features[8]
  vsrBrd2yz = features[9]

  vsrCs = features[10]
  vsrCsd2x = features[11]
  vsrCsd2yz = features[12]

  return params[0]*vcoul + params[1]*vcould2x + params[2]*vcould4x + params[3]*vcould4yz \
                         + params[4]*vsrPb + params[5]*vsrPbd2x + params[6]*vsrPbd2yz \
                         + params[7]*vsrBr + params[8]*vsrBrd2x + params[9]*vsrBrd2yz \
                         + params[10]*vsrCs + params[11]*vsrCsd2x + params[12]*vsrCsd2yz + params[13] + params[14]

def spsig_model(features,*params):
  x,y,z = features
  return params[0]*np.exp(params[1]*x + params[2]*(y**2 + z**2)) + params[3]

def ppsig_model(features,*params):
  x,y,z = features
  return params[0]*np.exp(params[1]*x + params[2]*(y**2 + z**2)) + params[3]

def pppi_model(features,*params):
  x,y,z = features
  return params[0]*np.exp(params[1]*x + params[2]*y**2 + params[3]*z**2) + params[4]

def ppn_br_to_pb_model(features,*params):
  x,y,z = features
  return params[0]*y*np.exp(params[1]*x + params[2]*z**2)

def ppn_pb_to_br_model(features,*params):
  x,y,z = features
  return params[0]*y*np.exp(params[1]*x + params[2]*z**2)

def spn_pb_to_br_model(features,*params):
  x,y,z = features
  return params[0]*y*np.exp(params[1]*x + params[2]*z**2)

def socpb_model(features,*params):
  vcoul = features[0]
  vcould2x = features[1]
  vcould2y = features[2]
  vcould4x = features[3]
  vcould4y = features[4]
  vcould4z = features[5]

  vsrPb = features[6]
  vsrPbd2x = features[7]
  vsrPbd2y = features[8]
  vsrPbd2z = features[9]

  vsrBr = features[10]
  vsrBrd2x = features[11]
  vsrBrd2y = features[12]
  vsrBrd2z = features[13]

  vsrCs = features[14]
  vsrCsd2x = features[15]
  vsrCsd2y = features[16]
  vsrCsd2z = features[17]

  return params[0]*vcoul + params[1]*vcould2x + params[2]*vcould2y \
                         + params[3]*vcould4x + params[4]*vcould4y + params[5]*vcould4z \
                         + params[6]*vsrPb + params[7]*vsrPbd2x + params[8]*vsrPbd2y +  params[9]*vsrPbd2z \
                         + params[10]*vsrBr + params[11]*vsrBrd2x + params[12]*vsrBrd2y + params[13]*vsrBrd2z\
                         + params[14]*vsrCs + params[15]*vsrCsd2x + params[16]*vsrCsd2y + params[17]*vsrCsd2z + params[18]

def socbr_par_model(features,*params):
  vcoul = features[0]
  vcould2x = features[1]
  vcould2y = features[2]
  vcould4x = features[3]
  vcould4y = features[4]
  vcould4z = features[5]

  vsrPb = features[6]
  vsrPbd2x = features[7]
  vsrPbd2y = features[8]
  vsrPbd2z = features[9]

  vsrBr = features[10]
  vsrBrd2x = features[11]
  vsrBrd2y = features[12]
  vsrBrd2z = features[13]

  vsrCs = features[14]
  vsrCsd2x = features[15]
  vsrCsd2y = features[16]
  vsrCsd2z = features[17]

  return params[0]*vcoul + params[1]*vcould2x + params[2]*vcould2y \
                         + params[3]*vcould4x + params[4]*vcould4y + params[5]*vcould4z \
                         + params[6]*vsrPb + params[7]*vsrPbd2x + params[8]*vsrPbd2y +  params[9]*vsrPbd2z \
                         + params[10]*vsrBr + params[11]*vsrBrd2x + params[12]*vsrBrd2y + params[13]*vsrBrd2z\
                         + params[14]*vsrCs + params[15]*vsrCsd2x + params[16]*vsrCsd2y + params[17]*vsrCsd2z + params[18]

def socbr_perp_model(features,*params):
  vcoul = features[0]
  vcould2x = features[1]
  vcould2y = features[2]
  vcould4x = features[3]
  vcould4y = features[4]
  vcould4z = features[5]

  vsrPb = features[6]
  vsrPbd2x = features[7]
  vsrPbd2y = features[8]
  vsrPbd2z = features[9]

  vsrBr = features[10]
  vsrBrd2x = features[11]
  vsrBrd2y = features[12]
  vsrBrd2z = features[13]

  vsrCs = features[14]
  vsrCsd2x = features[15]
  vsrCsd2y = features[16]
  vsrCsd2z = features[17]

  return params[0]*vcoul + params[1]*vcould2x + params[2]*vcould2y \
                         + params[3]*vcould4x + params[4]*vcould4y + params[5]*vcould4z\
                         + params[6]*vsrPb + params[7]*vsrPbd2x + params[8]*vsrPbd2y +  params[9]*vsrPbd2z \
                         + params[10]*vsrBr + params[11]*vsrBrd2x + params[12]*vsrBrd2y + params[13]*vsrBrd2z\
                         + params[14]*vsrCs + params[15]*vsrCsd2x + params[16]*vsrCsd2y + params[17]*vsrCsd2z + params[18]


def set_models(lut):
  lut.set("hoptype", "pbs", "model",  pbs_model)
  lut.set("hoptype", "pbp", "model",  pbp_model)
  lut.set("hoptype", "brp_par", "model",  brp_par_model)
  lut.set("hoptype", "brp_perp", "model",  brp_perp_model)
  
  lut.set("hoptype", "spsig", "model",  spsig_model)
  lut.set("hoptype", "ppsig", "model",  ppsig_model)
  lut.set("hoptype", "pppi", "model",  pppi_model)
  
  lut.set("hoptype", "ppn_br_to_pb", "model",  ppn_br_to_pb_model)
  lut.set("hoptype", "ppn_pb_to_br", "model",  ppn_pb_to_br_model)
  lut.set("hoptype", "spn_pb_to_br", "model",  spn_pb_to_br_model)
  
  lut.set("hoptype", "socpb", "model",  socpb_model)
  lut.set("hoptype", "socbr_perp", "model",  socbr_perp_model)
  lut.set("hoptype", "socbr_par", "model",  socbr_par_model)

