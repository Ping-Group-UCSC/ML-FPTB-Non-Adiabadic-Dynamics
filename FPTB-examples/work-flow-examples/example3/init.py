import sys
sys.path.insert(1,"../../common")
import gen_features_common as gfc

def makelut(atomindexdat, orbindexdat):
  lut = gfc.lookuptable()
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

    if orb_typ=="Br:p":
      ddir = ll[5]
      pdir = ll[6]
      spin = ll[7]
      lut.set("orb", orb_idx, "dir", ddir)
      lut.set("orb", orb_idx, "pdir", pdir)
      lut.set("orb", orb_idx, "spin", spin)
    
      if ddir=='x':
        atom_typ="Br1"
      if ddir=='y':
        atom_typ="Br2"
      if ddir=='z':
        atom_typ="Br3"
    
    lut.set("orb", orb_idx, "atom_typ", atom_typ)
    lut.set("orb", orb_idx, "ix", ix)
    lut.set("orb", orb_idx, "iy", iy)
    lut.set("orb", orb_idx, "iz", iz)

    atom_idx = lut.get("atom_typ",atom_typ,"ix",ix,"iy",iy,"iz",iz,"atom")[0]
    lut.set("orb", orb_idx, "atom", atom_idx)


  orb_index_file.close()
  return lut



