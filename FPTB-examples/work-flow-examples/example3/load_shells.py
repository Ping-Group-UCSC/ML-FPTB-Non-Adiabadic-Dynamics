import numpy as np


def set_shells(lut,NX,NY,NZ):
  mshells_max = 2
  
  pbpb_shells = np.loadtxt("../example3/shells/pbpb_shells.dat",dtype=int)
  pbbr_shells = np.loadtxt("../example3/shells/pbi_shells.dat",dtype=int)
  pbcs_shells = np.loadtxt("../example3/shells/pbma_shells.dat",dtype=int)
  br1pb_shells = np.loadtxt("../example3/shells/i1pb_shells.dat",dtype=int)
  br2pb_shells = np.loadtxt("../example3/shells/i2pb_shells.dat",dtype=int)
  br3pb_shells = np.loadtxt("../example3/shells/i3pb_shells.dat",dtype=int)
  br1cs_shells = np.loadtxt("../example3/shells/i1ma_shells.dat",dtype=int)
  br2cs_shells = np.loadtxt("../example3/shells/i2ma_shells.dat",dtype=int)
  br3cs_shells = np.loadtxt("../example3/shells/i3ma_shells.dat",dtype=int)
  br1br_shells = np.loadtxt("../example3/shells/i1i_shells.dat",dtype=int)
  br2br_shells = np.loadtxt("../example3/shells/i2i_shells.dat",dtype=int)
  br3br_shells = np.loadtxt("../example3/shells/i3i_shells.dat",dtype=int)
  brpb_shells = [br1pb_shells,br2pb_shells,br3pb_shells]
  brcs_shells = [br1cs_shells,br2cs_shells,br3cs_shells]
  brbr_shells = [br1br_shells,br2br_shells,br3br_shells]
  
  n_pbpb_shells = np.loadtxt("../example3/shells/n_pbpb_shells.dat",dtype=int)
  n_pbbr_shells = np.loadtxt("../example3/shells/n_pbi_shells.dat",dtype=int)
  n_pbcs_shells = np.loadtxt("../example3/shells/n_pbma_shells.dat",dtype=int)
  n_br1pb_shells = np.loadtxt("../example3/shells/n_i1pb_shells.dat",dtype=int)
  n_br2pb_shells = np.loadtxt("../example3/shells/n_i2pb_shells.dat",dtype=int)
  n_br3pb_shells = np.loadtxt("../example3/shells/n_i3pb_shells.dat",dtype=int)
  n_br1cs_shells = np.loadtxt("../example3/shells/n_i1ma_shells.dat",dtype=int)
  n_br2cs_shells = np.loadtxt("../example3/shells/n_i2ma_shells.dat",dtype=int)
  n_br3cs_shells = np.loadtxt("../example3/shells/n_i3ma_shells.dat",dtype=int)
  n_br1br_shells = np.loadtxt("../example3/shells/n_i1i_shells.dat",dtype=int)
  n_br2br_shells = np.loadtxt("../example3/shells/n_i2i_shells.dat",dtype=int)
  n_br3br_shells = np.loadtxt("../example3/shells/n_i3i_shells.dat",dtype=int)
  n_brpb_shells = [n_br1pb_shells,n_br2pb_shells,n_br3pb_shells]
  n_brcs_shells = [n_br1cs_shells,n_br2cs_shells,n_br3cs_shells]
  n_brbr_shells = [n_br1br_shells,n_br2br_shells,n_br3br_shells]
  
  
  for ix in range(NX):
    for iy in range(NY):
      for iz in range(NZ):
  
        #print ("%s %s %s" % (ix,iy,iz))
        atom1 = lut.get("atom_typ","Pb","ix",ix+1,"iy",iy+1,"iz",iz+1,"atom" )[0]
  
        #pb-pb shells
        for ishell in range(mshells_max): 
          j1shell = n_pbpb_shells[ishell+1] 
          j2shell = n_pbpb_shells[ishell+2] 
  
          for jshell in range(j1shell,j2shell):
            dx,jx = np.divmod( ix + pbpb_shells[jshell,0] , NX)
            dy,jy = np.divmod( iy + pbpb_shells[jshell,1] , NY)
            dz,jz = np.divmod( iz + pbpb_shells[jshell,2] , NZ)
            atom2 = lut.get("atom_typ","Pb","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]
  
            lut.set("atom", atom1, "pbpb_shell_%s" % (ishell+1), (atom2,dx,dy,dz))
  
        #pb-br shells
        for ishell in range(mshells_max): 
          j1shell = n_pbbr_shells[ishell+0] 
          j2shell = n_pbbr_shells[ishell+1] 
  
          for jshell in range(j1shell,j2shell):
            dx,jx = np.divmod( ix + pbbr_shells[jshell,0] , NX)
            dy,jy = np.divmod( iy + pbbr_shells[jshell,1] , NY)
            dz,jz = np.divmod( iz + pbbr_shells[jshell,2] , NZ)
            jtyp = pbbr_shells[jshell,3]
            atom2 = lut.get("atom_typ","Br%s" % jtyp,"ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]
  
            lut.set("atom", atom1, "pbbr_shell_%s" % (ishell+1), (atom2,dx,dy,dz))
  
        #pb-cs shells
        for ishell in range(mshells_max): 
          j1shell = n_pbcs_shells[ishell+0] 
          j2shell = n_pbcs_shells[ishell+1] 
  
          for jshell in range(j1shell,j2shell):
            dx,jx = np.divmod( ix + pbcs_shells[jshell,0] , NX)
            dy,jy = np.divmod( iy + pbcs_shells[jshell,1] , NY)
            dz,jz = np.divmod( iz + pbcs_shells[jshell,2] , NZ)
            atom2 = lut.get("atom_typ","Cs","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]
  
            lut.set("atom", atom1, "pbcs_shell_%s" % (ishell+1), (atom2,dx,dy,dz))
  
        
        for ityp in range(1,4):
          atom1 = lut.get("atom_typ","Br%s" % ityp,"ix",ix+1,"iy",iy+1,"iz",iz+1,"atom" )[0]
          #br-pb shells
          for ishell in range(mshells_max): 
            j1shell = n_brpb_shells[ityp-1][ishell+0] 
            j2shell = n_brpb_shells[ityp-1][ishell+1] 
  
            for jshell in range(j1shell,j2shell):
              dx,jx = np.divmod( ix + brpb_shells[ityp-1][jshell,0] , NX)
              dy,jy = np.divmod( iy + brpb_shells[ityp-1][jshell,1] , NY)
              dz,jz = np.divmod( iz + brpb_shells[ityp-1][jshell,2] , NZ)
              atom2 = lut.get("atom_typ","Pb","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]
  
              lut.set("atom", atom1, "brpb_shell_%s" % (ishell+1), (atom2,dx,dy,dz))
  
          #br-cs shells
          for ishell in range(mshells_max): 
            j1shell = n_brcs_shells[ityp-1][ishell+0] 
            j2shell = n_brcs_shells[ityp-1][ishell+1] 
  
            for jshell in range(j1shell,j2shell):
              dx,jx = np.divmod( ix + brcs_shells[ityp-1][jshell,0] , NX)
              dy,jy = np.divmod( iy + brcs_shells[ityp-1][jshell,1] , NY)
              dz,jz = np.divmod( iz + brcs_shells[ityp-1][jshell,2] , NZ)
              atom2 = lut.get("atom_typ","Cs","ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]
  
              lut.set("atom", atom1, "brcs_shell_%s" % (ishell+1), (atom2,dx,dy,dz))
  
          #br-br shells
          for ishell in range(mshells_max): 
            j1shell = n_brbr_shells[ityp-1][ishell+1] 
            j2shell = n_brbr_shells[ityp-1][ishell+2] 
  
            for jshell in range(j1shell,j2shell):
  
              dx,jx = np.divmod( ix + brbr_shells[ityp-1][jshell,0] , NX)
              dy,jy = np.divmod( iy + brbr_shells[ityp-1][jshell,1] , NY)
              dz,jz = np.divmod( iz + brbr_shells[ityp-1][jshell,2] , NZ)
              jtyp = brbr_shells[ityp-1][jshell,3]
              atom2 = lut.get("atom_typ","Br%s" % jtyp,"ix",jx+1,"iy",jy+1,"iz",jz+1,"atom" )[0]
  
              lut.set("atom", atom1, "brbr_shell_%s" % (ishell+1), (atom2,dx,dy,dz))


