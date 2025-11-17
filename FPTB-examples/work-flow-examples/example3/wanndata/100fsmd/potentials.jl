using LinearAlgebra 

qs = []

open("atom_index.dat") do f
  while ! eof(f)
    ll = readline(f)
    atom_type = split(ll)[2]
    if atom_type=="Cs"
      push!(qs,1.0)
    elseif atom_type=="Pb"
      push!(qs,2.0)
    elseif atom_type=="Br1" || atom_type=="Br2" || atom_type=="Br3"
      push!(qs,-1.0)
    end
  end
end

nq = size(qs,1)
LX = 12.00445150246208 #supercell dimensions
LY = 12.00445150246208
LZ = 12.00445150246208
cell = diagm([LX,LY,LZ])

ts = transpose(readdlm("atompos.dat"))
ts = cell * ts
        

#ewald summation of potential (not total energy)
# i = index of atom to calculate potential at
# see micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf
function ewald_phi(i, cell,qs,ts,sig,rcut,gcut):
  shortrange = 0.0
  longrange = 0.0
  selfint = 0.0

  gvecs = 2*pi*inv(cell)
  vol = abs(det(cell))
  nq = size(qs,1)

  for ix=-rcut:rcut
    for iy=-rcut:rcut
      for iz=-rcut:rcut
        for j=1:nq

          if(i==j && ix==0 && iy==0 && iz==0):
            continue
          end

          rdiff = norm(ts[:,i]-ts[:,j]+ cell*[ix,iy,iz])
          shortrange = shortrange + qs[j]/rdiff * erfc(rdiff/sqrt(2.0)/sig)
        end
      end
    end
  end

  for gx=-gcut:gcut
    for gy=-gcut:gcut
      for gz=-gcut:gcut

        if(gx==0 && gy==0 && gz==0)
          continue
        end

        gvec = vec(transpose([gx,gy,gz])*gvecs)
        for j=1:nq
          longrange = longrange + 2.0*pi/vol * exp(-0.5*sig^2 * norm(gvec)^2)/(norm(gvec)^2) * qs[j]
        end
      end
    end
  end


  selfint =  sqrt(2.0/pi)/sig * qs[i]

  return shortrange+longrange-selfint
end

for i=1:nq
  ewald_phi(i, cell,qs,ts,48,2,2)
end
