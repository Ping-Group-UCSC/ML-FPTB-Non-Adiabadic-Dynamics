
#writes a file which contains the orbital indices in the same order as Wannier90's orde

fout = open("orb_index.dat",'w')

nx = 4 
ny = 4 
nz = 4 

idx = 1

for typ in ["Pb:s", "Pb:p", "Br:p"]:
  for iz in range(nx):
    for iy in range(ny):
      for ix in range(nz):
        if typ == "Br:p":
          for d in ["z","y","x"]:
            for pdir in ["Z","X","Y"]:
              for spin in ["up","dn"]:
                fout.write("%s %s %s %s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1,d, pdir, spin) )
                idx = idx +1
        elif typ == "Pb:s":
          for spin in ["up","dn"]:
            fout.write("%s %s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1, spin) )
            idx = idx +1
        else:
          for pdir in ["Z","X","Y"]:
            for spin in ["up","dn"]:
              fout.write("%s %s %s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1, pdir, spin) )
              idx = idx +1




fout.close()
