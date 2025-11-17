
#writes a file which contains the atom indices in the same order as David's perovskite_structure.py : create_MAPbI3_cubic_cell

fout = open("atom_index.dat",'w')

nx = 2
ny = 2
nz = 2

idx = 1
for ix in range(nx):
  for iy in range(ny):
    for iz in range(nz):
      fout.write("%s Pb %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s I1 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s I2 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s I3 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s H_N1 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1
      
      fout.write("%s H_N2 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s H_N3 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s N %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s C %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s H_C1 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1
      
      fout.write("%s H_C2 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

      fout.write("%s H_C3 %s %s %s\n" % (idx, ix+1, iy+1, iz+1))
      idx = idx +1

fout.close()
