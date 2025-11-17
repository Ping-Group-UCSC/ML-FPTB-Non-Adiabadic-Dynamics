
#writes a file which contains the atom indices in the same order as David's perovskite_structure.py : create_MAPbI3_cubic_cell

fout = open("atom_index.dat",'w')

nx = 4 
ny = 4 
nz = 4 

idx = 1
for typ in ['Cs', 'Pb', 'Br1', 'Br2', 'Br3']:
  for iz in range(nx):
    for iy in range(ny):
      for ix in range(nz):
        fout.write("%s %s %s %s %s\n" % (idx, typ, ix+1, iy+1, iz+1))
        idx = idx +1

fout.close()
