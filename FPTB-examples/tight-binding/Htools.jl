module Htools

using Printf
using LinearAlgebra
using Arpack
using DelimitedFiles

export solvehamk, calcmmn, calcamn, projham

#find the bandstructure of a hamiltonian
#hammaker(k) = hamiltonian
#read the wavefunctions using read(wfnfile,Complex128,(dim(hamiltonian),nbands))
#gives evecs[ig,ibands]
# Options:
# savewfns: whether to save the wavefunctions
# wfnfmt : print in binary (false) or ascii (true)
# dir: output directory
# spars: use sparse matrix diagonalization, 
# nev: number of eigenvalues to find in sparse matrix diagonalization
# etarget: find eigenvalues close to this in sparse matrix diagonalization
#
#
# if(wfnfmt), the wfn file has the format
# re(band1,g1) im(band1,g1) re(band2,g1) im(band2,g1) ...
# re(band1,g2) im(band1,g2) re(band2,g2) im(band2,g2) ...
# ...
#
# Additionally, prints to file
# k.dat: kpoints in invang
# bands.dat: kpoints and energies
function solvehamk(hammaker, kpts, bands; savewfns=true,wfnfmt=false,verbose=true, dir="work", spars=false,nev=8,etarget=1.0,outfile="bands.dat")

  run(`mkdir -p $dir`)
  ndim = size(kpts,2)
  
  fk = open("$dir/k.dat","w")
  for ik = 1:size(kpts,1)
    kpt = vec(kpts[ik,:])
    for idim = 1:ndim
      @printf(fk,"%s ",kpt[idim])
    end
    @printf(fk,"\n")
  end
  close(fk)

  fout = open("$dir/$outfile","w")
  for ik = 1:size(kpts,1)
    if(verbose)
      println("doing ik: ",ik)
    end
    kpt = vec(kpts[ik,:])
    if(savewfns)
      run(`mkdir -p $dir/k$ik`)
    end
    if(spars)
      ens, evecs = eigs(hammaker(kpt); nev=nev, sigma=etarget)
      #ens = real(ens)
      bands = [1:nev]
    else
      ens, evecs = eigen(hammaker(kpt))
    end
    perm = sortperm(ens,by=real)
    if(savewfns)
      writedlm("$dir/k$ik/ens", ens[perm][bands], ' ' ) 
      evecs_sorted = evecs[:,perm]
      if(wfnfmt)
        nh = size(evecs_sorted,1)
        nb = size(bands,1)
        evecs_towrite = zeros(typeof(real(evecs_sorted[1,1])),(nh,2*nb))
        for ig=1:nh
          for ib=1:nb
            ib1 = bands[ib]
            evecs_towrite[ig,2*(ib-1)+1]=real(evecs_sorted[ig,ib1])
            evecs_towrite[ig,2*(ib-1)+2]=imag(evecs_sorted[ig,ib1])
          end
        end
        writedlm("$dir/k$ik/evecs", evecs_towrite, ' ' )
      else
        wfnfile = open("$dir/k$ik/evecs","w")
        write(wfnfile, evecs_sorted[:,bands])
        close(wfnfile)
      end
    end
    for b in sort(real(ens))
      for idim = 1:ndim
        @printf(fout,"%s ",kpt[idim])
      end
      @printf(fout,"%s\n",b)
    end
  end
  close(fout)
end

#project hamiltonian onto a subspace given by a projection operator
function projham(h,p,en)
  n = size(h,1)
  #u,s,v = svd(p)
  #tol = 1e-9
  #m = count(s,x->(abs(x)>tol)) # dim of subspace to project to

  q = eye(n) - p

  return p*h*p + p*h*q * inv(en*eye(n)-q*h*q) * q*h*p
end


#calculates overlap matrix elements
# ans[m,n,ik,ib] = < psi_k_m | psi_k+b_n >
# fw is a wavefunction transformer that tells
# what to do when BZ boundary is crossed.
# e.g. scrambleg for planewaves.
# fw(wfn,g) = wfn'.
#
# Input parameters:
# nbands = number of bands
# nkpts = number of kpoints
# nn = number of nearest neighbor kpoints
# nhs = dim of hilbert space
# nnkp = nearest neighbor kpoints data in
#        wanner90 format. Assumes that nnkp[i,1] is
#        increasing with i.
function calcmmn(bands, nbandstot,nkpts,nn,nhs,nnkp,fw; dir="work")
  nbands = size(bands,1)
  ans = zeros(ComplexF64,(nbands,nbands,nkpts,nn))
  for ik=1:nkpts
    wfn1file = open("$dir/k$ik/evecs","r")
    wfn1 = read(wfn1file,ComplexF64,(nhs,nbandstot))
    close(wfn1file)
    for ib=1:nn
      ikb = nnkp[(ik-1)*nn+ib,2]
      dg = vec(nnkp[(ik-1)*nn+ib,3:end])
      wfn2file = open("$dir/k$ikb/evecs","r")
      wfn2 = read(wfn2file,ComplexF64,(nhs,nbandstot))
      close(wfn2file)
      for m1=1:nbands
        for n1=1:nbands
          m = bands[m1]
          n = bands[n1]
          ans[m1,n1,ik,ib] = dot(wfn1[:,m], fw(wfn2[:,n],dg))
        end
      end
    end
  end
  return ans
end

#wannfns[:,n] = guess for the n-th wannier function
#ans[m,n,ik] = < psi_m_k | wann_n >
function calcamn(bands,nkpts,nwan,nhs,wannfns;dir="work")
  nbands = size(bands,1)
  ans = zeros(ComplexF64,(nbands,nwan,nkpts))
  for ik=1:nkpts
    wfnfile = open("$dir/k$ik/evecs","r")
    wfn = read(wfnfile,ComplexF64,(nhs,nbands))
    close(wfnfile)
    for m1=1:nbands
      for n=1:nwan
        m = bands[m1]
        ans[m,n,ik] = dot(wfn[:,m],wannfns[:,n])
      end
    end
  end
  return ans
end


end
