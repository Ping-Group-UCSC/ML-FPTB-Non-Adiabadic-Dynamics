import sys
sys.path.insert(1,"../../common")
import gen_features_common as gfc

import standard_positions as stdpos

import numpy as np

def set_phases(lut,NX,NY,NZ):

  for hoptype in ["pbs","pbp","brp_par","brp_perp","ppsig","pppi","ppn_br_to_pb", "ppn_pb_to_br" ,"spn_pb_to_br"] :
    for hopping in lut.get("hoptype",hoptype,"hopping"):
      lut.set("hopping",hopping,"phase",1.0)
  
  #spsig
  for hopping in lut.get("hoptype","spsig","hopping"):
    dx,dy,dz,orb1,orb2 = hopping
    atom1 = lut.get("orb",orb1,"atom")[0]
    atom2 = lut.get("orb",orb2,"atom")[0]
  
    ddir2 = lut.get("orb",orb2,"dir")[0]
    pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
    dr = stdpos.posvec0(lut,atom2)-stdpos.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
    if dr[0]>0.25 or dr[1]>0.25 or dr[2]>0.25:
      lut.set("hopping",hopping,"phase",1.0)
    else:
      lut.set("hopping",hopping,"phase",-1.0)
  
  
  
  #soc
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
  
          lut.set( "hopping",(0,0,0,orb1,orb3),"phase", 1.0j)
          lut.set( "hopping",(0,0,0,orb2,orb4),"phase", -1.0j )
          lut.set( "hopping",(0,0,0,orb2,orb5),"phase", -1.0  )
          lut.set( "hopping",(0,0,0,orb1,orb6),"phase", 1.0 )
          lut.set( "hopping",(0,0,0,orb3,orb6),"phase", 1.0j )
          lut.set( "hopping",(0,0,0,orb4,orb5),"phase", 1.0j )
  
        for atom2 in lut.get("atom_typ","Br1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
  
          orb1 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","X","orb")[0]
          orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","X","orb")[0]
          orb3 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Y","orb")[0]
          orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Y","orb")[0]
          orb5 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Z","orb")[0]
          orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Z","orb")[0]
  
          lut.set( "hopping",(0,0,0,orb1,orb3),"phase", 1.0j)
          lut.set( "hopping",(0,0,0,orb2,orb4),"phase", -1.0j )
          lut.set( "hopping",(0,0,0,orb2,orb5),"phase", -1.0  )
          lut.set( "hopping",(0,0,0,orb1,orb6),"phase", 1.0 )
          lut.set( "hopping",(0,0,0,orb3,orb6),"phase", 1.0j )
          lut.set( "hopping",(0,0,0,orb4,orb5),"phase", 1.0j )
  
  
        for atom2 in lut.get("atom_typ","Br2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
  
          orb1 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","X","orb")[0]
          orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","X","orb")[0]
          orb3 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Y","orb")[0]
          orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Y","orb")[0]
          orb5 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Z","orb")[0]
          orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Z","orb")[0]
  
          lut.set( "hopping",(0,0,0,orb1,orb3),"phase", 1.0j)
          lut.set( "hopping",(0,0,0,orb2,orb4),"phase", -1.0j )
          lut.set( "hopping",(0,0,0,orb2,orb5),"phase", -1.0  )
          lut.set( "hopping",(0,0,0,orb1,orb6),"phase", 1.0 )
          lut.set( "hopping",(0,0,0,orb3,orb6),"phase", 1.0j )
          lut.set( "hopping",(0,0,0,orb4,orb5),"phase", 1.0j )
  
  
        for atom2 in lut.get("atom_typ","Br3","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
  
          orb1 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","X","orb")[0]
          orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","X","orb")[0]
          orb3 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Y","orb")[0]
          orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Y","orb")[0]
          orb5 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Z","orb")[0]
          orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Z","orb")[0]
  
          lut.set( "hopping",(0,0,0,orb1,orb3),"phase", 1.0j)
          lut.set( "hopping",(0,0,0,orb2,orb4),"phase", -1.0j )
          lut.set( "hopping",(0,0,0,orb2,orb5),"phase", -1.0  )
          lut.set( "hopping",(0,0,0,orb1,orb6),"phase", 1.0 )
          lut.set( "hopping",(0,0,0,orb3,orb6),"phase", 1.0j )
          lut.set( "hopping",(0,0,0,orb4,orb5),"phase", 1.0j )

  

