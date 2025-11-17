
module Mymath

#include("mytools.jl")
using Mytools
using LinearAlgebra
using Combinatorics

export linspace, logspace
export almosteq, modrange, fermidirac, boseein, lor, mygauss, mylogit, mymod1, iround
export myslice
export trapz, trapzn, rieint
export projto, projto1, projaway, projaway1, gramschmidt, unitfy, myeig, matexp
export lerp, lerpcube
#export euex
export randgauss, randc, randcs, randc1, randcs1, coin,randrs,randrs1,fyshuf
export orbgen, sumperm

#linearly spaced array, endpoint inclusive by default
function linspace(start,stop,n;incl=true)
  if incl
    return collect(range(start,length=n,stop=stop))
  else
    return collect(range(start,length=n,stop=start+(stop-start)*(n-1)/n))
  end

end

#log spaced array, endpoint inclusive by default
# dense near the start of the array, sparse near the end
function logspace(start,stop,n;incl=true)
  return map(exp,linspace(log(start),log(stop),n,incl=incl))
end



function almosteq(d::Float64,x::Float64,y::Float64)
  #approximate equality
  return abs(x-y)<d
end

#modulus for reals
#returns something in the range [lb,ub)
function modrange(lb,ub,x)
  return mod(x-lb,ub-lb) + lb
end

function lor(x,c,b)
  return (1.0/pi)*b/((x-c)*(x-c)+b*b)
end

#gaussian, normalized
function mygauss(x,mu,sig)
  return (1.0/sig/sqrt(2.0*pi))*exp(-0.5*((x-mu)/sig)^2)
end

#logistic function
function mylogit(x,mu,k)
  return 1.0/(1.0+exp(-1.0*k*(x-mu)))
end

function fermidirac(x,mu,kt)
  return 1.0/(1.0+exp((x-mu)/kt))
end

function boseein(x,mu,kt)
  return 1.0/(exp((x-mu)/kt)-1.0)
end



#
# Integers
#

#double factorial
function factorial2(n)
  ans = n
  m = n
  while true
    m = m-2
    if m>0
      ans = ans*m
    else
      break
    end
  end
  return ans
end


#1-based modulus
function mymod1(x,m)
  return mod((x-1),m) + 1
end

function iround(x)
  return convert(Int,round(x))
end

#
# Linear algebra utilities
#

#project v onto u
function projto(u::Array{T,1},v::Array{T,1}) where T
  return dot(u,v)/dot(u,u) * u
end

#project v onto u but do not return a vector
function projto1(u::Array{T,1},v::Array{T,1}) where T
  return dot(u,v)/dot(u,u) 
end

#component of v perp to u
function projaway(u::Array{T,1},v::Array{T,1}) where T
  return v-projto(u,v)
end

#component of v perp to u but do not return a vector
function projaway1(u::Array{T,1},v::Array{T,1}) where T
  return norm(projaway(u,v))
end


#normalize a vector
function unitfy(v)
  return (1.0/norm(v)) * v
end

#input: vectors vs[i,:]
#output: orthonormalized vectors ans[i,:]
function gramschmidt(vs::Array{T,2}) where T
  n = size(vs,1)
  m = size(vs,2)
  ans = zeros(typeof(vs[1,1]),(n,m))
  for i=1:n
    temp = vec(vs[i,:])
    for j=1:i-1
      temp = temp - projto(vec(ans[j,:]),vec(vs[i,:]))
    end
    nn = norm(temp)
    if(nn>0)
      ans[i,:] = temp/nn
    else
      ans[i,:] = 0
    end
  end
  return ans
end

#gram-schmidt with error threshold for 
#removing linearly dependent vectors
function gramschmidt(vs::Array{T,2},th) where T
  n = size(vs,1)
  m = size(vs,2)
  ans = zeros(typeof(vs[1,1]),(0,m))
  for i=1:n
    temp = vec(vs[i,:])
    for j=1:size(ans,1)
      temp = temp - projto(vec(ans[j,:]),vec(vs[i,:]))
    end
    nn = norm(temp)
    if(nn>th)
      ans = cat(ans,transpose(temp/nn),dims=1)
    end
  end
  return ans
end



#project a vector onto a subspace
#spanned by vectors ss[i,:]
function projto(ss::Array{T,2}, v::Array{T,1}) where T
  if size(ss,1)==1
    return projto(vec(ss),v)
  else
    ss1 = gramschmidt(ss)
    ans = 0 * v
    for i=1:size(ss,1)
      ans = ans + projto(ss1[i,:],v)
    end
    return ans
  end
end


#project a vector away from a subspace
#spanned by vectors ss[i,:]
function projaway(ss::Array{T,2}, v::Array{T,1}) where T
  return v-projto(ss,v)
end


#diagonalize a matrix, sorting eigenvalues and eigenvectors in 
#increasing order of Re(eigenvalues).
#Eigenvectors are normalized: evecs[:,i] is the ith eigenvector
function myeig(m)
  ens, evecs = eig(m)
  perm = sortperm(ens,by=real)
  return (ens[perm],evecs[:,perm])
end

#matrix exponential
function matexp(m::Array{T,2}) where T
  ens, evecs = eig(m)
  eens = zeros(typeof(ens[1]),size(m))
  for ii=1:size(ens,1)
    eens[ii,ii] = exp(ens[ii])
  end
  return evecs * eens * inv(evecs)
end


#
# Numerical integration
#

#trapezoidal rule numerical integration
function trapz(xs,ys)
  n = min(size(xs,1), size(ys,1))
  ans = 0
  for i = 1:(n-1)
    ans = ans+(xs[i+1]-xs[i])*(ys[i+1]+ys[i])/2.0
  end
  return ans
end

#riemann integral
function rieint(xs,ys)
  n = min(size(xs,1), size(ys,1))
  ans = 0
  for i = 1:(n-1)
    ans = ans+(xs[i+1]-xs[i])*ys[i]
  end
  return ans
end

# array utility for slicing n dim arrays
# returns n-m dim array.
# myslice(arr,i1,i2) = arr[:,:, ... , i1, i2]
function myslice(x,ii...)
  n = ndims(x)
  m = length(ii)
  y = gensym()
  defy = eval(:($y = $x))
  aa = cat( [y], [:(:) for i=1:(n-m)], [ii...], dims=1)
  return eval(Expr(:ref, aa...))
end


#multidimensional trapezoidal integration
# on a rectangular (but not regular) grid
# hs[ix,iy] = value of intrgrand at xs[ix],ys[iy]
function trapzn(hs,xs...)
  nd = length(xs) 
  if(nd ==1)
    return trapz(xs[1],hs)
  else
    n = size(xs[nd],1)
    newxs = xs[1:nd-1]
    temph = [trapzn(myslice(hs,i),newxs...) for i=1:n]
    return trapz(xs[nd],temph)
  end
end

#
# Interpolation
#

#linear interpolation
# xs and ys will be automatically sorted according to xs
function lerp(xs0,ys0,xp)

  perm = sortperm(xs0, by=real)
  xs = xs0[perm]
  ys = ys0[perm]

  j = binsearchl(xs,xp)

  if(j==0)
    if(xs[1]>xp)
      return ys[1]
    else
      return ys[end]
    end
  end

  x0 = xs[j]
  x1 = xs[j+1]
  y0 = ys[j]
  y1 = ys[j+1]
  
  return y0 + (y1-y0)*(xp-x0)/(x1-x0)
end

#interpolate a point within a n-dimensional hypercube
#coordinates [0,1]
function lerpcube(cube,p)
  ndim = size(p,1)
  if(ndim==1)
    return cube[1]*(1.0-p[1])+cube[2]*p[1]
  else
    ars = fill(2,ndim-1)
    cube1 = zeros(Float64,tuple(ars...))
    mf = multifor(ars)
    for j=1:2^(ndim-1)
      idxs = vec(mf[j,:]) .+ 1
      x = p[end]
      y0 = geti(cube,cat(idxs,[1],dims=1))
      y1 = geti(cube,cat(idxs,[2],dims=1))
      val = y0*(1-x)+y1*x
      seti!(cube1,idxs,val)
    end
    return lerpcube(cube1,p[1:end-1])
  end
end



#
# Differential equations
#


#
# Random numbers
#

#gaussian random distribution,
#using Box-muller transform
function randgauss(mu,sig)
  u1 = rand()
  u2 = rand()
  return mu+ sig*(sqrt(-2.0*log(u1))*cos(2*pi*u2));
end

# generate random variable on [0,1] from
# inverse cumulative distribution function F^-1
function rand01(finv)
  fx = rand()
  return finv(fx)
end

# random complex phase
function randc1()
  return exp(1.0im*rand()*2*pi)
end

#random complex number less than norm 1
function randc()
  return rand01(sqrt)*randc1()
end

#random complex vector with norm1
function randcs1(n)
  z = [randc() for i=1:n]
  return z/norm(z)
end


#random complex n-vector,
#with norm < 1
function randcs(n)
  return rand01(u->u^(1/2/n))*randcs1(n) 
end


#uniform random real n-vector on n-ball
function randrs(n)
  while true
    v = [rand()*2.0-1.0 for i=1:n]
    if norm(v)<1
      return v
      break
    end
  end
end

#uniform random real n-vector on n-1 sphere
function randrs1(n)
  v = randrs(n)
  return v/norm(v)
end


#coin flip with some probability
function coin(p)
  if (rand()<p)
    return true
  else
    return false
  end
end

#Fisher-Yates shuffle, 
#generates random permutation of [1:n]
function fyshuf(n)
  a = [1:n]
  for i =1:n-1
    j = rand(i:n)
    temp = a[i]
    a[i] = a[j]
    a[j] = temp
  end
  return a
end

#
# Group Theory
#

#"plain vanilla" orbit generation
# gen[i](w) = w'
# weq(w,w') = (is w==w' ?)
function orbgen(w,gens,weq)
  orb = [w]
  que = [w]
  while(size(que,1)>0)
    wp = pop!(que)
    for g in gens
      wpp = g(wp)

      norb = size(orb,1)
      found = false
      for ix=1:norb
        found = weq(orb[ix],wpp)
        if found
          break
        end
      end

      if !(found)
        push!(orb,wpp)
        push!(que,wpp)
      end
    end
  end

  return orb        
end

#return \sum_P f(P(xs)) 
# where the sum is over all permutations of the argument list xs
# the "plus" operator in the summation can be arbitrary
function sumperm(f,xs,plusop)
  ps = permutations(xs)
  next = iterate(ps)
  (i0,state0) = next
  accu = f(i0...)
  next = iterate(ps,state0)
  while next!==nothing
    (i,state) = next
    accu = plusop(accu, f(i...))
    next = iterate(ps,state)
  end
  return accu
end


end
