import numpy as np


def get_vsr(i, cell, qs, ts, msk, la, rcut):
  '''
    i:     atom number
    cell:  cell dimentions in 3,3 array
    qs:    atomic charges
    ts:    atomic positions
    msk:   atom type mask
    la:    crystal radii of ion
    rcut:  cutoff for periodic images

    Returns:
        vsr:  short range potential
  '''

  vsr = 0.0
  vol = abs(np.linalg.det(cell))
  nq = len(qs)
  charge = qs[msk][0] 
  for ix in range(-rcut,rcut+1):
    for iy in range(-rcut,rcut+1):
      for iz in range(-rcut,rcut+1):
        for j in range(nq):

          if(i==j and ix==0 and iy==0 and iz==0):
            continue
       
          if not msk[j]:
            continue
          charge = qs[j]
          rdiff = np.linalg.norm(ts[:,i]-ts[:,j]+ np.matmul(cell,np.array([ix,iy,iz])))
          vsr += (1./(8.*np.pi*(la**3)))*np.exp(-1.0*rdiff/la)

  return vsr


def d2vsr(rdiffx,rdiffy,rdiffz,la):

  rdiff = np.linalg.norm([rdiffx, rdiffy, rdiffz])

  return np.exp(-rdiff/la)*(rdiffx**2/(8.*la**5*rdiff**2) - (rdiffy**2+rdiffz**2)/(8.*la**4*rdiff**(3./2.)))


def get_d2vsr(i, cell, qs, ts, msk, la, rcut):
  '''
    i:     atom number
    cell:  cell dimentions in 3,3 array
    qs:    atomic charges
    ts:    atomic positions
    msk:   atom type mask
    la:    crystal radii of ion
    rcut:  cutoff for periodic images

    Returns:
        vsrd2x, vsrd2y vsrd2z:  second derivative of short range potential for x, y, z
  '''

  vsr2dx = 0.0
  vsr2dy = 0.0
  vsr2dz = 0.0
  vol = abs(np.linalg.det(cell))
  nq = len(qs)
  charge = qs[msk][0] 
  for ix in range(-rcut,rcut+1):
    for iy in range(-rcut,rcut+1):
      for iz in range(-rcut,rcut+1):
        for j in range(nq):

          if(i==j and ix==0 and iy==0 and iz==0):
            continue
       
          if not msk[j]:
            continue
          charge = qs[j]
          rdiffx, rdiffy, rdiffz = np.absolute(ts[:,i]-ts[:,j] + np.matmul(cell,np.array([ix,iy,iz])))
          #if rdiffx > 0.:
          vsr2dx += d2vsr(rdiffx, rdiffy, rdiffz, la)
          #if rdiffy > 0.:
          vsr2dy += d2vsr(rdiffy, rdiffz, rdiffx, la)
          #if rdiffz > 0.:
          vsr2dz += d2vsr(rdiffz, rdiffx, rdiffy, la)

  return np.array([vsr2dx, vsr2dy, vsr2dz])

def d2coul(rdiffx, rdiffy, rdiffz, q):
   
  rdiff = np.linalg.norm([rdiffx, rdiffy, rdiffz])

  return -q * (2.*rdiffx**2 + rdiffy**2 + rdiffz**2)/(4.0*np.pi*rdiff**(5.0/2.0))

def get_d2coul(i, cell, qs, ts, rcut):
  '''
    i:     atom number
    cell:  cell dimentions in 3,3 array
    qs:    atomic charges
    ts:    atomic positions
    rcut:  cutoff for periodic images

    Returns:
        vcould2x, vcould2y vcould2z:  second derivative of coulomb potential for x, y, z
  '''
  epsil = 1
  vcoul2dx = 0.0
  vcoul2dy = 0.0
  vcoul2dz = 0.0
  vol = abs(np.linalg.det(cell))
  nq = len(qs)
  for ix in range(-rcut,rcut+1):
    for iy in range(-rcut,rcut+1):
      for iz in range(-rcut,rcut+1):
        for j in range(nq):

          if(i==j and ix==0 and iy==0 and iz==0):
            continue
       
          charge = qs[j]
          rdiffx, rdiffy, rdiffz = np.absolute(ts[:,i]-ts[:,j] + np.matmul(cell,np.array([ix,iy,iz])))
          #if rdiffx > 0.:
          vcoul2dx += d2coul(rdiffx, rdiffy, rdiffz, charge)
          #if rdiffy > 0.:
          vcoul2dy += d2coul(rdiffy, rdiffz, rdiffx, charge)
          #if rdiffz > 0.:
          vcoul2dz += d2coul(rdiffz, rdiffx, rdiffy, charge)

  return np.array([vcoul2dx, vcoul2dy, vcoul2dz])


def d4coul(rdiffx, rdiffy, rdiffz, q):
   
  rdiff = np.linalg.norm([rdiffx, rdiffy, rdiffz])

  return 3.0 * q * (8.*rdiffx**4 - (24.0*rdiffx**2*(rdiffy**2 + rdiffz**2)) + 3.0*(rdiffy**2 + rdiffz**2)**2)/(4.0*np.pi*rdiff**(9.0/2.0))

def get_d4coul(i, cell, qs, ts, rcut):
  '''
    i:     atom number
    cell:  cell dimentions in 3,3 array
    qs:    atomic charges
    ts:    atomic positions
    rcut:  cutoff for periodic images

    Returns:
        vcould4x, vcould4y vcould4z: fourth derivative of coulomb potential for x, y, z
  '''
  epsil = 1
  vcoul4dx = 0.0
  vcoul4dy = 0.0
  vcoul4dz = 0.0
  vol = abs(np.linalg.det(cell))
  nq = len(qs)
  for ix in range(-rcut,rcut+1):
    for iy in range(-rcut,rcut+1):
      for iz in range(-rcut,rcut+1):
        for j in range(nq):

          if(i==j and ix==0 and iy==0 and iz==0):
            continue
       
          charge = qs[j]
          rdiffx, rdiffy, rdiffz = np.absolute(ts[:,i]-ts[:,j] + np.matmul(cell,np.array([ix,iy,iz])))
          #if rdiffx > 0.:
          vcoul4dx += d4coul(rdiffx, rdiffy, rdiffz, charge)
          #if rdiffy > 0.:
          vcoul4dy += d4coul(rdiffy, rdiffz, rdiffx, charge)
          #if rdiffz > 0.:
          vcoul4dz += d4coul(rdiffz, rdiffx, rdiffy, charge)

  return np.array([vcoul4dx, vcoul4dy, vcoul4dz])




