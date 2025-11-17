
#writes a file which contains the orbital indices in the same order as David's order

fout = open("orb_index.dat",'w')

nx = 2
ny = 2
nz = 2

idx = 1
for ix in range(nx):
  for iy in range(ny):
    for iz in range(nz):
      for typ in ["Pb:s", "Pb:p", "I:p"]:
        if typ == "Pb:s":
          for spin in ["up","dn"]:
            fout.write("%s %s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1, spin) )
            idx = idx +1

        elif typ == "Pb:p":
          for pdir in ["X","Y","Z"]:
            for spin in ["up","dn"]:
              fout.write("%s %s %s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1, pdir, spin) )
              idx = idx +1

        elif typ == "I:p":
          for d in ["x","y","z"]:
            for pdir in ["X","Y","Z"]:
              for spin in ["up","dn"]:
                fout.write("%s %s %s %s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1,d, pdir, spin) )
                idx = idx +1

fout.close()
