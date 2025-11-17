# This script re-introduces the signs of spsig and nnii_cage
# which have been removed for the ML training step.
# Run this before running the tight-binding calculation
#
#

import sys
import numpy as np

sys.path.append("../gen_features/")
import gen_features_common as gfc

# load the hoppings to be fixed
hopfilename = "hopping.dat.new"
hops = np.loadtxt(hopfilename)

NX = 2
NY = 2
NZ = 2

#
# initialize lookuptables
#
lut =  gfc.makelut("../gen_features/atom_index.dat","../gen_features/orb_index.dat") 

fixedhops = []

for hh in hops:
    dx,dy,dz,orb1,orb2 = hh[0:5].astype(int)

    print(lut.get("orb",orb1,"orb_typ")[0], lut.get("orb",orb2,"orb_typ")[0] )

    #sp
    if lut.get("orb",orb1,"orb_typ")[0] == "Pb:s" and lut.get("orb",orb2,"orb_typ")[0] == "I:p" :
        atom1 = lut.get("orb",orb1,"atom")[0]
        atom2 = lut.get("orb",orb2,"atom")[0]

        ddir2 = lut.get("orb",orb2,"dir")[0]
        pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
        if ddir2==pdir2: #this is spsig
            dr = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
            spsign = -1.0
            if dr[0]>0.25 or dr[1]>0.25 or dr[2]>0.25:
                spsign = 1.0
            fixedhops.append([dx,dy,dz,orb1,orb2, spsign*hh[5], 0.0 ])
        else: #this is spn_pb_to_i : signs don't need to be fixed
            fixedhops.append([dx,dy,dz,orb1,orb2, hh[5], hh[6] ])

    #I p - I p 
    elif lut.get("orb",orb1,"orb_typ")[0] == "I:p" and lut.get("orb",orb2,"orb_typ")[0] == "I:p" :

        atom1 = lut.get("orb",orb1,"atom")[0]
        atom2 = lut.get("orb",orb2,"atom")[0]

        ddir1 = lut.get("orb",orb1,"dir")[0]
        pdir1 = lut.get("orb",orb1,"pdir")[0].lower()

        ddir2 = lut.get("orb",orb2,"dir")[0]
        pdir2 = lut.get("orb",orb2,"pdir")[0].lower()
        #this is nnii_cage
        if ddir1==pdir2 and ddir2==pdir1 and ddir1!=pdir1 and ddir2!=pdir2 and atom1!=atom2:
            dr = gfc.posvec0(lut,atom2)-gfc.posvec0(lut,atom1)+np.array([NX*dx,NY*dy,NZ*dz])
            nnii_cage_sign = -1.0
            if ddir1=="x":
                idx1 = 0
            elif ddir1=="y":
                idx1 =1
            elif ddir1=="z":
                idx1 =2
            if ddir2=="x":
                idx2 = 0
            elif ddir2=="y":
                idx2 =1
            elif ddir2=="z":
                idx2 =2

            if dr[idx1]*dr[idx2]>0:
                nnii_cage_sign = 1.0
            fixedhops.append([dx,dy,dz,orb1,orb2, nnii_cage_sign*hh[5], 0.0 ])
        elif pdir1!=pdir2 and atom1==atom2: #this is I SOC
            big = abs(hh[5])
            small = hh[6]

            if lut.get("orb",orb1,"spin")[0] =="up" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="up" and lut.get("orb",orb2,"pdir")[0] =="Y":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, -1.0*big ])

            if lut.get("orb",orb1,"spin")[0] =="dn" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="dn" and lut.get("orb",orb2,"pdir")[0] =="Y":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, 1.0*big ])

            if lut.get("orb",orb1,"spin")[0] =="dn" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="up" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, -1.0*big, small ])

            if lut.get("orb",orb1,"spin")[0] =="up" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="dn" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, 1.0*big, small ])

            if lut.get("orb",orb1,"spin")[0] =="up" and lut.get("orb",orb1,"pdir")[0] =="Y" and lut.get("orb",orb2,"spin")[0] =="dn" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, -1.0*big ])

            if lut.get("orb",orb1,"spin")[0] =="dn" and lut.get("orb",orb1,"pdir")[0] =="Y" and lut.get("orb",orb2,"spin")[0] =="up" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, -1.0*big ]) 


            
        else:
            fixedhops.append([dx,dy,dz,orb1,orb2, hh[5], hh[6] ])
                
    #Pb p - Pb p 
    elif lut.get("orb",orb1,"orb_typ")[0] == "Pb:p" and lut.get("orb",orb2,"orb_typ")[0] == "Pb:p" :

        atom1 = lut.get("orb",orb1,"atom")[0]
        atom2 = lut.get("orb",orb2,"atom")[0]

        pdir1 = lut.get("orb",orb1,"pdir")[0].lower()
        pdir2 = lut.get("orb",orb2,"pdir")[0].lower()

        if pdir1!=pdir2 and atom1==atom2: #this is Pb SOC
            big = abs(hh[5])
            small = hh[6]

            if lut.get("orb",orb1,"spin")[0] =="up" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="up" and lut.get("orb",orb2,"pdir")[0] =="Y":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, -1.0*big ])

            if lut.get("orb",orb1,"spin")[0] =="dn" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="dn" and lut.get("orb",orb2,"pdir")[0] =="Y":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, 1.0*big ])

            if lut.get("orb",orb1,"spin")[0] =="dn" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="up" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, -1.0*big, small ])

            if lut.get("orb",orb1,"spin")[0] =="up" and lut.get("orb",orb1,"pdir")[0] =="X" and lut.get("orb",orb2,"spin")[0] =="dn" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, 1.0*big, small ])

            if lut.get("orb",orb1,"spin")[0] =="up" and lut.get("orb",orb1,"pdir")[0] =="Y" and lut.get("orb",orb2,"spin")[0] =="dn" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, -1.0*big ])

            if lut.get("orb",orb1,"spin")[0] =="dn" and lut.get("orb",orb1,"pdir")[0] =="Y" and lut.get("orb",orb2,"spin")[0] =="up" and lut.get("orb",orb2,"pdir")[0] =="Z":
                fixedhops.append([dx,dy,dz,orb1,orb2, small, -1.0*big ]) 

            
        else:
            fixedhops.append([dx,dy,dz,orb1,orb2, hh[5], hh[6] ])
         

    else:
        fixedhops.append([dx,dy,dz,orb1,orb2, hh[5], hh[6] ])

np.savetxt(hopfilename+"_fixed",np.array(fixedhops), fmt=("%4i","%4i","%4i","%4i","%4i","%8.4f","%8.4f"))
