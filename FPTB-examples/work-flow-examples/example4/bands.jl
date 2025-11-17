using DelimitedFiles
using TheoryTools.Htools
using TheoryTools.TB
using TheoryTools.Kspace

hoppingdat = readdlm("hopping.dat")
#hoppingdat = readdlm("orig_hopping.dat")

kpath0 = [ 0.0 0.0 0.0 ; 0.5 0.0 0.0 ; 0.5 0.5 0.0]
kpath = kpts_path(kpath0,10)

solvehamk(k->tbhammaker3d(k,208,hoppingdat), kpath, [], savewfns=false, verbose=true, spars= false)
