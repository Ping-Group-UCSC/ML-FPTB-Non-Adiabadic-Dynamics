module Mytb

#include("mytools.jl")
using Mytools
using Mymath
using LinearAlgebra
using SparseArrays

export   tbhammaker1d, tbhammaker2d, tbhammaker3d, tb1d


#
# The central datastructure in this tight binding code is hopsdat,
# an array of hoppings between pairs of orbitals.
# Each row of hopsdat = [h,j,l, from, to, re(val), im(val)],
# where the hopping is between orbitals "from" and "to", which is in unit cells displaced by (h,j,l).
# h,j,l = 000 means orbitals "from" and "to" are in the same unit cell.
# h,j,l = 100 means orbitals "from" and "to" are in different unit cells, 1 lattice constant apart in first lattice direction.
# re(val), im(val) are the real and imaginary parts of the hopping.
#
# Once created, hopsdat allows construction of the tight-binding Hamiltonian at any k point, using the tbhammaker functions.

function tbhammaker1d(kpt1, nhs, hopsdat)
  #kpt1 = [kpt]
  #kpt in crystal coords (in units of recip latt. vecs
  #no need to specify inverse hoppings
  ans = zeros(ComplexF64,(nhs,nhs))
  kpt = kpt1[1]

  for ih in 1:size(hopsdat,1)
    dcell = iround(hopsdat[ih,1])
    from = iround(hopsdat[ih,2])
    to = iround(hopsdat[ih,3])
    val = hopsdat[ih,4] + 1.0im*hopsdat[ih,5]

    phase = exp(2.0im *pi* dcell*kpt)
    #println(kpt)
    ans[from,to] = ans[from,to] + phase*val
    if(dcell!=0 || from!=to)
      ans[to,from] = ans[to,from] + conj(phase*val)
    end

  end
  return ans
end

#atomic positions are specified to fix the phases
#atpos[i] = position of atom i in crystal coords
#kpt1 = [kpt]
#kpt in crystal coords (in units of recip latt. vecs
#no need to specify inverse hoppings
function tb1d(kpt1,atpos,hopsdat)
  nhs = size(atpos,1)
  ans = zeros(ComplexF64,(nhs,nhs))
  kpt = kpt1[1]

  for ih in 1:size(hopsdat,1)
    dcell = iround(hopsdat[ih,1])
    from = iround(hopsdat[ih,2])
    to = iround(hopsdat[ih,3])
    val = hopsdat[ih,4] + 1.0im*hopsdat[ih,5]

    phase = exp(2.0im *pi* (dcell + atpos[to] - atpos[from]) *kpt)
    #println(kpt)
    ans[from,to] = ans[from,to] + phase*val
    if(dcell!=0 || from!=to)
      ans[to,from] = ans[to,from] + conj(phase*val)
    end
  end
  return ans
end



function tbhammaker2d(kpt, nhs, hopsdat)
  #kpt in crystal coords (in units of recip latt. vecs
  ans = zeros(ComplexF64,(nhs,nhs))

  for ih in 1:size(hopsdat,1)
    dcell = [iround(hopsdat[ih,i]) for i = 1:2]
    from = iround(hopsdat[ih,3])
    to = iround(hopsdat[ih,4])
    val = hopsdat[ih,5] + 1.0im*hopsdat[ih,6]

    phase = exp(2.0im *pi* dot(dcell,kpt))
    #println(from,", ",to)
    ans[from,to] = ans[from,to] + phase*val
    if(dcell[1]!=0 || dcell[2]!=0 || from!=to)
      ans[to,from] = ans[to,from] + conj(phase*val)
    end

  end
  return ans
end

function tbhammaker3d(kpt, nhs, hopsdat; spars=false)
  #kpt in crystal coords (in units of recip latt. vecs
  #spas = constuct a sparse matrix

  if(spars)
    nhops = size(hopsdat,1)

    rcidxs = zeros(Int,nhops)
    vals = zeros(ComplexF64,nhops)
    rcidxs1 = Int[]
    vals1 = Complex{Float64}[]

    for ih in 1:nhops
      dcell = [round(Int64,hopsdat[ih,i]) for i = 1:3]
      from = round(Int64,hopsdat[ih,4])
      to = round(Int64,hopsdat[ih,5])
      val = hopsdat[ih,6] + 1.0im*hopsdat[ih,7]

      phase = exp(2.0im *pi* dot(dcell,kpt))
      #println(from,", ",to)

      rcidxs[ih] = flati1([nhs,nhs],[from-1,to-1])
      vals[ih] = phase*val
      if(dcell[1]!=0 || dcell[2]!=0 || dcell[3]!=0 || from!=to)

        push!(rcidxs1, flati1([nhs,nhs],[to-1,from-1]))
        push!(vals1,  conj(phase*val))
      end

    end

    rcidxs2 = cat(rcidxs,rcidxs1,dims=1)
    vals2 = cat(vals,vals1,dims=1)
    nhops2 = size(vals2,1)

    perm = sortperm(rcidxs2)
    rcidxss = rcidxs2[perm]
    valss = vals2[perm]

    rcidxu = zeros(Int,nhops2)
    rowidxu = zeros(Int,nhops2)
    colidxu = zeros(Int,nhops2)
    valsu = zeros(ComplexF64,nhops2)
    jj=1
    currnum = rcidxss[1]
    rcidxu[1] = currnum
    valsu[1] = valss[1]
    for ii=2:nhops2
      if(currnum<rcidxss[ii])
        currnum = rcidxss[ii]
        jj=jj+1
        rcidxu[jj] = currnum
      end
      valsu[jj] = valsu[jj] + valss[ii]
    end

    for ii=1:jj
      rc = multiforhelp1([nhs,nhs],rcidxu[ii])
      rowidxu[ii] = rc[1]+1
      colidxu[ii] = rc[2]+1
    end

    #println(sort(rowidxu[1:jj]))
    #println(colidxu[1:jj])
    #println(valsu[1:jj])

    ans = sparse(rowidxu[1:jj],colidxu[1:jj],valsu[1:jj])
  else

    ans = zeros(ComplexF64,(nhs,nhs))

    for ih in 1:size(hopsdat,1)
      dcell = [round(Int64,hopsdat[ih,i]) for i = 1:3]
      from = round(Int64,hopsdat[ih,4])
      to = round(Int64,hopsdat[ih,5])
      to = iround(hopsdat[ih,5])
      val = hopsdat[ih,6] + 1.0im*hopsdat[ih,7]

      phase = exp(2.0im *pi* dot(dcell,kpt))
      #println(from,", ",to)
      ans[from,to] = ans[from,to] + phase*val
      if(dcell[1]!=0 || dcell[2]!=0 || dcell[3]!=0 || from!=to)
        ans[to,from] = ans[to,from] + conj(phase*val)
      end

    end
  end
  return ans
end

end
