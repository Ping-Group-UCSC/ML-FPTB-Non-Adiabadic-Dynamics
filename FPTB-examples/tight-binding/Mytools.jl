module Mytools


export multiforhelp, multii_c0, multii_c1, multiforhelp1, multifor, flati1, flati, flati_c1, flati_c0, flati_f0, flati_f1, geti, seti!
export flatten,map2, aoa2d, nub,binsearch, binsearchl
export indicvec, deleten


#
#some functions related to indexing
#

# b = "bases"
# n = ans[bl] 
#   + ans[bl-1]*b[bl]
#   + ans[bl-2]*b[bl]*b[bl-1]
#   +  ...
#zero-indexed
# It's the inverse of flati.
function multiforhelp(b,n)
  bl = size(b,1)
  ans = zeros(Int,bl)
  temp = n
  for i=1:bl
    q = div(temp,b[bl-i+1])
    r = mod(temp,b[bl-i+1])
    temp = q
    ans[bl-i+1] = r
  end
  return ans
end

#alias for multiforhelp: cstyle ordering, and 0-indexed
multii_c0 = multiforhelp 

#1-indexed
function multii_c1(b,n)
  return multii_c0(b, n-1) .+ 1
end

# This is like multifor_help, 
# except first index varies fastest.
# It's the inverse of flati1.
# n = ans[1] 
#   + ans[2]*b[1] 
#   + ans[3]*b[1]*b[2] + ...
#zero-indexed
function multiforhelp1(b,n)
  bl = size(b,1)
  ans = zeros(Int,bl)
  temp = n
  for i=1:bl
    q = div(temp,b[i])
    r = mod(temp,b[i])
    temp = q
    ans[i] = r
  end
  return ans
end

# such that 0<= ans[j,i] < b[i] 
# last index varies fastest
#zero-indexed
function multifor(b)
  bl = size(b,1)
  anslen = int(prod(b))
  ans = zeros(Int,(anslen,bl))
  for i = 0:anslen-1
    ans[i+1,:] = multiforhelp(b,i)
  end
  return ans
end



#collapses multi-index to single index
#fortran-style ordering: first index varies fastest
#  for example, for a 2D matrix:
#  (i_1,i_2) -> i_1 + i_2 * n_1 
#zero-indexed
function flati1(arrshape,arri)
  n = size(arri,1)
  ans = 0
  for ii=1:n
    temp = arri[ii]
    for jj=2:ii
      temp = temp * arrshape[jj-1] 
    end
    ans = ans + temp
  end
  return ans
end

#collapses multi-index to single index
#c-style ordering: first index varies slowest
#  for example, for a 2D matrix:
#  (i_1,i_2) -> i_2 + i_1 * n_2 
#zero-indexed
function flati(arrshape,arri)
  return flati1(reverse(arrshape),reverse(arri))
end

#1-indexed version of flati
function flati_c1(arrshape,arri)
  return flati(arrshape, arri .- 1) +1
end

flati_c0 = flati
flati_f0 = flati1

#1-indexed version of flati1
function flati_f1(arrshape,arri)
  return flati1(arrshape, arri .- 1) +1
end




#geti(arr,[i1,i2,...]) = arr[i1,i2,...]
function geti(arr,vi)
  return arr[flati_f1(size(arr),vi)]
end

#seti!(arr,[i1,i2,...], val) 
#  means
#  arr[i1,i2,...] = val
function seti!(arr, arri, val)
  arr[flati_f1(size(arr), arri1)] = val
end

#
# Array utils
#

#expand the first level of an array of arrays
function flatten(xss)
  ans = []
  for xs in xss
    ans = cat(ans,xs,dims=1)
  end
  return ans
end

#map over 2 arrays
function map2(f,xs,ys)
  n = min(size(xs,1),size(ys,1))
  return [f(xs[i],ys[i]) for i=1:n]
end

#convert array of arrays into a 2d array
function aoa2d(aoa)
  n1 = size(aoa,1)
  n2 = size(aoa[1],1)
  ans = Array{typeof(aoa[1][1])}(undef,(n1,n2))
  for i=1:n1
    for j = 1:n2
      ans[i,j] = aoa[i][j]
    end
  end
  return ans
end

#tensor product of  multidimensional arrays
#function tensorprd(t1,t2)

#find unique elements of an array 
#answer in increasing order
function nub(arr)
  n = size(arr,1)
  t = typeof(arr[1])
  y = Array{t}(undef,n)
  x = sort(arr)

  jj=1
  currnum = x[1]
  y[1] = currnum
  for ii=2:n
    if(currnum<x[ii]) 
      currnum = x[ii]
      jj = jj+1
      y[jj] = currnum
    end
  end

  return y[1:jj]
end

#binary search of a sorted array
#binsearch(a,j) = i if a(i) = j
#binsearch(a,j) = 0 if j is not found
function binsearch(a,j)
  lena = size(a,1)
  lowlim = 1
  highlim = lena

  ans =0

  for ii=1:lena
    if(highlim<lowlim+2)
      if(a[lowlim] == j)
        ans = lowlim
      elseif(a[highlim]==j)
        ans = highlim
      else
        ans = 0
      end
      break
    else
      mid = div(lowlim+highlim,2)
      if(a[mid]<j) 
        lowlim=mid
      else
        highlim=mid
      end
    end
  end

  return ans
end

#binary search of a sorted array
#binsearchl(a,j) = i if a[i] <= j, a[i+1]>j 
#binsearchl(a,j) = 0 if a[1]>j or a[end]<j
function binsearchl(a,j)
  lena = size(a,1)
  lowlim = 1
  highlim = lena

  ans =0

  if(a[lowlim]>j || a[highlim]<j)
    return 0
  end

  for ii=1:lena
    if(highlim<lowlim+2)
      if(a[lowlim]<=j && a[highlim]>j)
        ans = lowlim
      else
        ans = 0
      end
      break
    else
      mid = div(lowlim+highlim,2)
      if(a[mid]<j) 
        lowlim=mid
      else
        highlim=mid
      end
    end
  end

  return ans
end

#create an "indicator vector"
#ans[i] = 0 if i!=j
#ans[j] = x
function indicvec(n,j,x)
  ans = zeros(typeof(x),n)
  ans[j] = x
  return ans
end

#removes the nth element of an array, 1-indexed
function deleten(v,n::Int)
  boolarr = fill(true,size(v,1))
  boolarr[n] = false
  return getindex(v,boolarr)
end

#removes some elements of an array
# ns = list of indexes to remove
function deleten(v,ns::Array{Int,1})
  boolarr = fill(true,size(v,1))
  for n in ns
    boolarr[n] = false
  end
  return getindex(v,boolarr)
end



# create a poor man's hash table
# the keys are vectors and the values are integers
# input: rr = 2D integer array
# output: pmh[i1,i2,...] = ii if rr[ii,:] = i1,i2,...
#   else, pmh[i1,i2,...] = 0
#function poormanhash(arr)

end
