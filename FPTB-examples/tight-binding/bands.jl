using Htools
using Mytb
using Kspace
using DelimitedFiles



# Import hopping data
hoppingdat = readdlm("hopping.dat")

#Specify cell dimensions here
nx = 2
ny = 2
nz = 2
ncells = nx*ny*nz
nh = ncells*26

# This makes a bandstructure along a k-path
kpath0 = [ 0.0 0.0 0.0 ; 0.5 0.5 0.5 ; 0.5 0.0 0.0 ; 0.0 0.0 0.0 ]
kpath = kpts_path(kpath0,20)


#use spars = true for sparse matrix diagonalization of nev bands, close to energy etarget
solvehamk(k->tbhammaker3d(k,nh,hoppingdat), kpath, [], savewfns=false,verbose=true, spars=false, nev=8,etarget=4.85 )



  


