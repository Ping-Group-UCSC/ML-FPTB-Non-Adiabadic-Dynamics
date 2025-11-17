import numpy as np

#
# construct hoppings of each type here
# dx, dy, dz, orb1, orb2

#use this to set only spin up hoppings and soc hoppings
def set_hoppings(lut,NX,NY,NZ):

  for i in lut.get("orb_typ","Pb:p","spin","up","orb"):
    lut.set("hoptype", "pbp", "hopping",(0,0,0,i,i))
  
  for i in lut.get("orb_typ","Pb:s","spin","up","orb"):
    lut.set("hoptype", "pbs", "hopping",(0,0,0,i,i))
  
  for orb1 in lut.get("orb_typ","Br:p","spin","up","orb"):
    ddir = lut.get("orb",orb1,"dir")[0]
    pdir = lut.get("orb",orb1,"pdir")[0].lower()
    if ddir==pdir:
      lut.set("hoptype", "brp_par", "hopping",(0,0,0,orb1,orb1))
    else:
      lut.set("hoptype", "brp_perp", "hopping",(0,0,0,orb1,orb1))
  
  for ix in range(NX):
    for iy in range(NY):
      for iz in range(NZ):
  
        atom1 = lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
        for neighbor in lut.get("atom", atom1, "pbbr_shell_1"): 
          atom2,dx,dy,dz = neighbor
  
          for orb1 in lut.get("atom",atom1,"spin","up","orb"):
            for orb2 in lut.get("atom",atom2,"spin","up","orb"):
  
              ddir2 = lut.get("orb",orb2,"dir")[0]
              pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
              if lut.get("orb",orb1,"orb_typ")[0] == "Pb:s":
                if ddir2==pdir2:
                  lut.set("hoptype", "spsig", "hopping",(dx,dy,dz,orb1,orb2))
                else:
                  lut.set("hoptype", "spn_pb_to_br", "hopping",(dx,dy,dz,orb1,orb2))
              elif lut.get("orb",orb1,"orb_typ")[0] == "Pb:p":
                pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
                if ddir2==pdir1:
                  if ddir2==pdir2:
                    lut.set("hoptype", "ppsig", "hopping",(dx,dy,dz,orb1,orb2))
                  else:
                    lut.set("hoptype", "ppn_pb_to_br", "hopping",(dx,dy,dz,orb1,orb2))
                else:
                  if pdir1==pdir2:
                    lut.set("hoptype", "pppi", "hopping",(dx,dy,dz,orb1,orb2))
                  elif ddir2==pdir2:
                    lut.set("hoptype", "ppn_br_to_pb", "hopping",(-dx,-dy,-dz,orb2,orb1))
  
  
  # sigma bonds involving p orbitals only
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
                  lut.set("hoptype", "nnpbpb", "hopping",(dx,dy,dz,orb1,orb2))
  
  
  #includes br-br nearest neighbor bonds involving p-orbitals perpendicular to Br-Pb 
  for ix in range(NX):
    for iy in range(NY):
      for iz in range(NZ):
  
        atom1 = lut.get("atom_typ","Br1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
        for neighbor in lut.get("atom", atom1, "brbr_shell_1"): 
          atom2,dx,dy,dz = neighbor
  
          for orb1 in lut.get("atom",atom1,"spin","up","orb"):
            for orb2 in lut.get("atom",atom2,"spin","up","orb"):
              pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
              pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  
              ddir1 = lut.get("orb",orb1,"dir")[0]
              ddir2 = lut.get("orb",orb2,"dir")[0]
              if pdir1!=pdir2 and pdir1==ddir2 and pdir2==ddir1:
                lut.set("hoptype", "nnbrbr_cage", "hopping",(dx,dy,dz,orb1,orb2))
  
        atom1 = lut.get("atom_typ","Br2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
        for neighbor in lut.get("atom", atom1, "brbr_shell_1"): 
          atom2,dx,dy,dz = neighbor
          atom_typ2 = lut.get("atom",atom2,"atom_typ")[0]
          if atom_typ2 == "Br3":
  
            for orb1 in lut.get("atom",atom1,"spin","up","orb"):
              for orb2 in lut.get("atom",atom2,"spin","up","orb"):
                pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
                pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
  
                ddir1 = lut.get("orb",orb1,"dir")[0]
                ddir2 = lut.get("orb",orb2,"dir")[0]
                if pdir1!=pdir2 and pdir1==ddir2 and pdir2==ddir1:
                  lut.set("hoptype", "nnbrbr_cage", "hopping",(dx,dy,dz,orb1,orb2))
  
  
  #sigma bonding "across" Br-Pb-Br only
  for ix in range(NX):
    for iy in range(NY):
      for iz in range(NZ):
  
        atom1 = lut.get("atom_typ","Br1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
        atom2 = lut.get("atom_typ","Br2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
        atom3 = lut.get("atom_typ","Br3","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom")[0]
  
        for atom4 in [atom1,atom2,atom3]:
  
          for neighbor in lut.get("atom", atom4, "brbr_shell_2"): 
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
  
              for orb4 in lut.get("atom",atom4,"spin","up","orb_typ","Br:p","orb"):
                for orb5 in lut.get("atom",atom5,"spin","up","orb_typ","Br:p","orb"):
                  pdir4 = lut.get("orb",orb4,"pdir")[0].lower()
                  pdir5 = lut.get("orb",orb5,"pdir")[0].lower()
  
                  ddir4 = lut.get("orb",orb4,"dir")[0]
                  ddir5 = lut.get("orb",orb5,"dir")[0]
                  if pdir4==pdir5 and pdir4==bdir and pdir4==ddir4 :
                    lut.set("hoptype", "nnbrbr_sig", "hopping",(dx,dy,dz,orb4,orb5))
  
  
  
  #socbr_par  #this is when one of the orbitals is parallel to the bond 
  #socbr_perp  #this is when both of the orbitals are perpendicular to the bond
  
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
  
          lut.set("hoptype", "socpb", "hopping",(0,0,0,orb1,orb3))
          lut.set("hoptype", "socpb", "hopping",(0,0,0,orb2,orb4))
          lut.set("hoptype", "socpb", "hopping",(0,0,0,orb2,orb5))
          lut.set("hoptype", "socpb", "hopping",(0,0,0,orb1,orb6))
          lut.set("hoptype", "socpb", "hopping",(0,0,0,orb3,orb6))
          lut.set("hoptype", "socpb", "hopping",(0,0,0,orb4,orb5))
  
        for atom2 in lut.get("atom_typ","Br1","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
  
          orb1 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","X","orb")[0]
          orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","X","orb")[0]
          orb3 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Y","orb")[0]
          orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Y","orb")[0]
          orb5 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Z","orb")[0]
          orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Z","orb")[0]
  
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb1,orb3))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb2,orb4))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb2,orb5))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb1,orb6))
          lut.set("hoptype", "socbr_perp", "hopping",(0,0,0,orb3,orb6))
          lut.set("hoptype", "socbr_perp", "hopping",(0,0,0,orb4,orb5))
  
  
        for atom2 in lut.get("atom_typ","Br2","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
  
          orb1 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","X","orb")[0]
          orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","X","orb")[0]
          orb3 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Y","orb")[0]
          orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Y","orb")[0]
          orb5 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Z","orb")[0]
          orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Z","orb")[0]
  
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb1,orb3))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb2,orb4))
          lut.set("hoptype", "socbr_perp", "hopping",(0,0,0,orb2,orb5))
          lut.set("hoptype", "socbr_perp", "hopping",(0,0,0,orb1,orb6))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb3,orb6))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb4,orb5))
  
  
        for atom2 in lut.get("atom_typ","Br3","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom"):
  
          orb1 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","X","orb")[0]
          orb2 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","X","orb")[0]
          orb3 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Y","orb")[0]
          orb4 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Y","orb")[0]
          orb5 = lut.get("atom",atom2,"spin","up","orb_typ","Br:p","pdir","Z","orb")[0]
          orb6 = lut.get("atom",atom2,"spin","dn","orb_typ","Br:p","pdir","Z","orb")[0]
  
          lut.set("hoptype", "socbr_perp", "hopping",(0,0,0,orb1,orb3))
          lut.set("hoptype", "socbr_perp", "hopping",(0,0,0,orb2,orb4))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb2,orb5))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb1,orb6))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb3,orb6))
          lut.set("hoptype", "socbr_par", "hopping",(0,0,0,orb4,orb5))


