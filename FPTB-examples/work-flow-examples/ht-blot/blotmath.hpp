#ifndef BLOTMATH
#define BLOTMATH

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <slepceps.h>
#include <Eigen/Dense>
#include <math.h>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include "blot/basics.hpp"

using namespace std::complex_literals;

//round-towards-minus-infinity integer division (aka knuth division, aka python remainder)
struct knuthdiv_t{
  int quot;
  int rem;
};

knuthdiv_t knuthdiv(int a, int b){ //works only for b!=0
  knuthdiv_t ans;
  auto dv = std::div(a,b);
  int sgnb = (b > 0) - (b < 0);
  int sgnrem = (dv.rem > 0 ) - (dv.rem < 0);
  ans.rem = dv.rem * sgnrem * sgnb;
  ans.quot = dv.quot - (sgnrem == -sgnb);
  return ans;
}


//
// Operations on std::vector
//

//scalar product of a std::vector
template<class T>
std::vector<T> operator* (std::vector<T> v, T s){
  std::vector<T> ans;
  for(int i=0; i<v.size(); ++i){
    ans.push_back(v[i]*s);
  }
  return ans;
}

template<class T>
std::vector<T> operator* (T s, std::vector<T> v){
  return v * s;
}

template<class T>
std::vector<T> operator/ (std::vector<T> v, T s){
  return v * (1/s);
}



//element wise product of std::vector
template<class T>
std::vector<T> operator* (std::vector<T> v1, std::vector<T> v2){
  std::vector<T> ans;
  for(int i=0; i<v1.size(); ++i){
    ans.push_back(v1[i]*v2[i]);
  }
  return ans;
}

//element wise sum of std::vector
template<class T>
std::vector<T> operator+ (std::vector<T> v1, std::vector<T> v2){
  std::vector<T> ans;
  for(int i=0; i<v1.size(); ++i){
    ans.push_back(v1[i]+v2[i]);
  }
  return ans;
}

template<class T>
std::vector<T> operator- (std::vector<T> v1, std::vector<T> v2){
  std::vector<T> ans;
  for(int i=0; i<v1.size(); ++i){
    ans.push_back(v1[i]-v2[i]);
  }
  return ans;
}



//broadcasted sum of std::vector
template<class T>
std::vector<T> operator+ (std::vector<T> v, T s){
  std::vector<T> ans;
  for(int i=0; i<v.size(); ++i){
    ans.push_back(v[i]+s);
  }
  return ans;
}

template<class T>
std::vector<T> operator+ (T s, std::vector<T> v){
  return v + s;
}


//
//
// PETSC WRAPPERS
//
//

PetscErrorCode PetscInitWrap(){
  PetscErrorCode ierr;

  int argc = 1;
  char** argv;
  argv = new char*[argc];
  argv[0] = (char*)"init wrapper";
  char help[] = "see blot/blotmath.cpp";
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);

  return ierr;
}


PetscErrorCode SlepcInitWrap(){
  PetscErrorCode ierr;

  int argc = 1;
  char** argv;
  argv = new char*[argc];
  argv[0] = (char*)"init wrapper";
  char help[] = "see blot/blotmath.cpp";
  ierr = SlepcInitialize(&argc,&argv,(char*)0,help); if(ierr) return ierr; 

  return ierr;
}



//wrappers around petsc mat and vec
struct MatWrap{
  int m;
  int n;
  Mat A;
  MatWrap(int m1, int n1);
  ~MatWrap();
  std::vector<PetscScalar> get();
  std::vector<PetscScalar> get(std::vector<int> s1, std::vector<int> s2);
  std::vector<PetscScalar> get(allseq as1, std::vector<int> s2);
  std::vector<PetscScalar> get(std::vector<int> s1, allseq as2);
  PetscScalar get(int t1, int t2);
};

MatWrap::MatWrap(int m1, int n1){
  m = m1;
  n = n1;
  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD,&A); //CHKERRQ(ierr);
  //ierr = MatCreate(PETSC_COMM_SELF,&A); //CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m1,n1);//CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);//CHKERRQ(ierr);
  ierr = MatSetUp(A);//CHKERRQ(ierr);
}

MatWrap::~MatWrap(){
  PetscErrorCode ierr;
  ierr = MatDestroy(&A);//CHKERRQ(ierr);
}

std::vector<PetscScalar> MatWrap::get(){
  PetscErrorCode ierr;
  std::vector<PetscScalar> ans(m*n, 0.0);
  auto s1 = idxseq(m);
  auto s2 = idxseq(n);
  int * a1 = &s1[0];
  int * a2 = &s2[0];
  ierr = MatGetValues(A,m,a1,n,a2, &ans[0]);
  return ans;
}

std::vector<PetscScalar> MatWrap::get(std::vector<int> s1, std::vector<int> s2){
  PetscErrorCode ierr;
  std::vector<PetscScalar> ans(s1.size()*s2.size(), 0.0);
  int * a1 = &s1[0];
  int * a2 = &s2[0];
  ierr = MatGetValues(A,m,a1,n,a2, &ans[0]);
  return ans;
}

std::vector<PetscScalar> MatWrap::get( allseq as1, std::vector<int> s2){
  PetscErrorCode ierr;
  std::vector<PetscScalar> ans(m*s2.size(), 0.0);
  auto s1 = idxseq(m);
  int * a1 = &s1[0];
  int * a2 = &s2[0];
  ierr = MatGetValues(A,m,a1,s2.size(),a2,&ans[0]);
  return ans;
}

std::vector<PetscScalar> MatWrap::get( std::vector<int> s1, allseq as2){
  PetscErrorCode ierr;
  std::vector<PetscScalar> ans(s1.size()*n, 0.0);
  auto s2 = idxseq(n);
  int * a1 = &s1[0];
  int * a2 = &s2[0];
  ierr = MatGetValues(A,s1.size(),a1,n,a2,&ans[0]);
  return ans;
}


PetscScalar MatWrap::get(int t1, int t2){
  PetscScalar v;
  MatGetValue(A,t1,t2,&v);
  return v;
}




struct VecWrap{
  int n;
  Vec V; 
  VecWrap(int n1);
  ~VecWrap();
  PetscErrorCode set(std::vector<PetscScalar> values);
  PetscErrorCode set(std::vector<int> s, std::vector<PetscScalar> values);
  PetscErrorCode set(int n, PetscScalar value);
  PetscScalar get(int n);
};

VecWrap::VecWrap(int n1){
    n = n1;
    PetscErrorCode ierr;
    ierr = VecCreate(PETSC_COMM_WORLD,&V); //CHKERRQ(ierr);
    //ierr = VecCreate(PETSC_COMM_SELF,&V); //CHKERRQ(ierr);
    ierr = VecSetSizes(V,PETSC_DECIDE,n1);//CHKERRQ(ierr);
    ierr = VecSetFromOptions(V);//CHKERRQ(ierr);
}


VecWrap::~VecWrap(){
  PetscErrorCode ierr;
  ierr = VecDestroy(&V);//CHKERRQ(ierr);
}


PetscErrorCode VecWrap::set(std::vector<PetscScalar> values){
  PetscErrorCode ierr;
  auto s1 = idxseq(n);
  int * a1 = &s1[0];
  VecSetValues(V,n,a1,&values[0],INSERT_VALUES);
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  return ierr;
}

PetscErrorCode VecWrap::set(std::vector<int> s, std::vector<PetscScalar> values){
  PetscErrorCode ierr;
  int * a1 = &s[0];
  VecSetValues(V,s.size(), a1 , &values[0],INSERT_VALUES);
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  return ierr;
}

PetscErrorCode VecWrap::set(int n, PetscScalar value){
  PetscErrorCode ierr;
  VecSetValue(V,n,value,INSERT_VALUES);
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  return ierr;
}

PetscScalar VecWrap::get(int n){
  PetscScalar val;
  const int n1 = n;
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
  VecGetValues(V, 1, &n, &val);
  return val;
}

//
// Linear Algebra
//

//petsc matrix inverse
PetscErrorCode invertWrap( MatWrap & A, MatWrap & Ainv){
  PetscErrorCode ierr;

  //copy the values into a block diagonal matrix 
  std::vector<PetscScalar> A_values = A.get();

  Mat A_baij;
  PetscInt idxm[1]={0}, idxn[1]={0};
  const PetscScalar * A_inv_values;
  MatCreateSeqBAIJ(PETSC_COMM_SELF,A.m,A.m,A.n,1,NULL,&A_baij);
  MatSetOption(A_baij,MAT_ROW_ORIENTED,PETSC_FALSE);
  MatSetValuesBlocked(A_baij,1,idxm,1,idxn,&A_values[0],INSERT_VALUES);
  MatAssemblyBegin(A_baij,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_baij,MAT_FINAL_ASSEMBLY);
  MatInvertBlockDiagonal(A_baij,&A_inv_values);

  //set values of Ainv  //check column major order?
  auto s1 = idxseq(Ainv.m);
  auto s2 = idxseq(Ainv.n);
  int * a1 = &s1[0];
  int * a2 = &s2[0];
  ierr = MatSetValues(Ainv.A,Ainv.m,a1,Ainv.n,a2,A_inv_values,INSERT_VALUES);
  MatAssemblyBegin(Ainv.A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Ainv.A,MAT_FINAL_ASSEMBLY);

  MatDestroy(&A_baij);


  return ierr;
}



//petsc multidim fft 
//ndim 	- the ndim-dimensional transform
//dim 	- array of size ndim, dim[i] contains the vector length in the i-dimension 

PetscErrorCode fftWrap(PetscInt ndim,const PetscInt dim[],PetscScalar * input,  PetscScalar * output){
  PetscErrorCode ierr;
  Mat A;
  Vec x,y,z;

  PetscInt ntot=1;
  for(auto ii=0;ii<ndim;++ii){
    ntot = ntot * dim[ii];
  }

  auto s1 = idxseq(ntot);
  int * a1 = &s1[0];

  MatCreateFFT(PETSC_COMM_WORLD,ndim,dim,MATFFTW,&A);
  MatCreateVecsFFTW(A,&x,&y,&z);
  VecSetValues(x,ntot,a1,input,INSERT_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  MatMult(A,x,y);
  
  VecGetValues(y,ntot,a1,output);

  VecDestroy(&x);
  VecDestroy(&y);
  VecDestroy(&z);
  MatDestroy(&A);

  return ierr;
}


//
// Ewald sums
//


//Ewald summation of total energy
//cell(all,i) is the ith lattice vector.
//qs = charges
//ts[i] = positions of charge i
//sig = erf parameter
//rcut, gcut = realspace/recipspace cutoffs

double ewald_tote(const Eigen::Matrix3d & cell, const std::vector<Eigen::Vector3d> & ts, const std::vector<double> & qs, double sig, int rcut, int gcut){

  double shortrange = 0.0;
  double longrange = 0.0;
  double selfint = 0.0;

  auto gvecs = 2* M_PI *cell.inverse();
  auto vol = abs(cell.determinant());

  int nq = qs.size();

  //int nt = ts.size();
  //int ndim1 = ts[0].size();
  //std::vector<Eigen::Matrix<double,1,Eigen::Dynamic>> tsvec;
  //for (auto ii=0; ii<nt; ++ii){
  //  Eigen::Matrix<double,1,Eigen::Dynamic> temp(1,ndim1);
  //  
  //  for (auto jj=0; jj<ndim1; ++jj){
  //    temp(jj) = ts[ii][jj];
  //  }
  //  tsvec.push_back(temp);
  //}

  for (auto ix=-rcut; ix<rcut+1; ++ix){
    for (auto iy=-rcut; iy<rcut+1; ++iy){
      for (auto iz=-rcut; iz<rcut+1; ++iz){
        for (auto ii=0; ii<nq; ++ii){
          for (auto jj=0; jj<nq; ++jj){

            if(ii==jj && ix==0 && iy==0 && iz==0) continue;

            Eigen::Vector3d dd {1.0*ix, 1.0*iy, 1.0*iz};

            auto rr = ts[ii] - ts[jj] + cell*dd;
            auto rdiff  = rr.norm();
            shortrange = shortrange + 0.5*qs[ii]*qs[jj]/rdiff * erfc(rdiff/sqrt(2.0)/sig);
          }
        }
      }
    }
  }

  for (auto gx=-gcut; gx<gcut+1; ++gx){
    for (auto gy=-gcut; gy<gcut+1; ++gy){
      for (auto gz=-gcut; gz<gcut+1; ++gz){

        if(gx==0 && gy==0 && gz==0) continue;
            
        Eigen::Matrix<double,1,3> gg {1.0*gx, 1.0*gy, 1.0*gz};

        auto gvec = gg * gvecs ;
        auto gdiff2 = gvec.dot(gvec);

        std::complex<double> strfac = 0.0;

        for (auto ii=0; ii<nq; ++ii){ 
          double tdotg = gvec * ts[ii];
          strfac = strfac+qs[ii]*exp(1.0i* tdotg );
        }

        longrange = longrange + 2.0*M_PI/vol * exp(-0.5*sig*sig * gdiff2)/gdiff2 * abs(strfac)*abs(strfac);
      }
    }
  }


  for (auto ii=0; ii<nq; ++ii){
    selfint = selfint + 1.0/sqrt(2.0*M_PI)/sig * qs[ii]*qs[ii];
  }

  return shortrange+longrange-selfint;
}

//ewald summation of potential (not total energy)
//ia = index of atom to calculate potential at
// see micro.stanford.edu/mediawiki/images/4/46/Ewald_notes.pdf
double ewald_phi(int ia, const Eigen::Matrix3d & cell, const std::vector<Eigen::Vector3d> & ts, const std::vector<double> & qs, double sig, int rcut, int gcut){

  double shortrange = 0.0;
  double longrange = 0.0;
  double selfint = 0.0;

  auto gvecs = 2* M_PI *cell.inverse();
  auto vol = abs(cell.determinant());

  int nq = qs.size();

  //int nt = ts.size();
  //int ndim1 = ts[0].size();
  //std::vector<Eigen::Matrix<double,1,Eigen::Dynamic>> tsvec;
  //for (auto ii=0; ii<nt; ++ii){
  //  Eigen::Matrix<double,1,Eigen::Dynamic> temp(1,ndim1);
  //  
  //  for (auto jj=0; jj<ndim1; ++jj){
  //    temp(jj) = ts[ii][jj];
  //  }
  //  tsvec.push_back(temp);
  //}

  for (auto ix=-rcut; ix<rcut+1; ++ix){
    for (auto iy=-rcut; iy<rcut+1; ++iy){
      for (auto iz=-rcut; iz<rcut+1; ++iz){
        for (auto jj=0; jj<nq; ++jj){

          if(ia==jj && ix==0 && iy==0 && iz==0) continue;

          Eigen::Matrix<double,3,1> dd {1.0*ix, 1.0*iy, 1.0*iz};

          auto rr = ts[ia] - ts[jj] + cell*dd;
          auto rdiff  = rr.norm();
          shortrange = shortrange + qs[jj]/rdiff * erfc(rdiff/sqrt(2.0)/sig);
        }
      }
    }
  }

  for (auto gx=-gcut; gx<gcut+1; ++gx){
    for (auto gy=-gcut; gy<gcut+1; ++gy){
      for (auto gz=-gcut; gz<gcut+1; ++gz){

        if(gx==0 && gy==0 && gz==0) continue;
            
        Eigen::Matrix<double,1,3> gg {1.0*gx, 1.0*gy, 1.0*gz};

        auto gvec = gg * gvecs;
        auto gdiff2 = gvec.dot(gvec);

        for (auto jj=0; jj<nq; ++jj){
          auto gr = gvec.dot(ts[ia]-ts[jj]);
          longrange = longrange + 4.0*M_PI/vol * cos(gr) * exp(-0.5*sig*sig * gdiff2)/gdiff2 * qs[jj];
        }
      }
    }
  }


  selfint = sqrt(2.0/M_PI)/sig * qs[ia];

  //cout << "ewald" << shortrange << " " << longrange << " " << selfint <<endl;

  return shortrange+longrange-selfint;
}


//ewald summation of force
//ia = index of atom to calculate efield at
// see http://www.pages.drexel.edu/~cfa22/msim/node49.html
Eigen::Vector3d ewald_field(int ia, const Eigen::Matrix3d & cell, const std::vector<Eigen::Vector3d> & ts, const std::vector<double> & qs, double sig, int rcut, int gcut){

  Eigen::Vector3d shortrange {0.0, 0.0, 0.0};
  Eigen::Vector3d longrange {0.0, 0.0, 0.0};

  auto gvecs = 2* M_PI *cell.inverse();
  auto vol = abs(cell.determinant());

  int nq = qs.size();

  for (auto ix=-rcut; ix<rcut+1; ++ix){
    for (auto iy=-rcut; iy<rcut+1; ++iy){
      for (auto iz=-rcut; iz<rcut+1; ++iz){
        for (auto jj=0; jj<nq; ++jj){

          if(ia==jj && ix==0 && iy==0 && iz==0) continue;

          Eigen::Matrix<double,3,1> dd {1.0*ix, 1.0*iy, 1.0*iz};

          auto rr = ts[ia] - ts[jj] + cell*dd;
          auto rdiff  = rr.norm();
          shortrange = shortrange + qs[jj]*rr/rdiff/rdiff * (sqrt(2/M_PI)/sig*exp(-0.5*rdiff*rdiff/sig/sig) +  (1.0/rdiff) * erfc(rdiff/sqrt(2.0)/sig));
        }
      }
    }
  }

  for (auto gx=-gcut; gx<gcut+1; ++gx){
    for (auto gy=-gcut; gy<gcut+1; ++gy){
      for (auto gz=-gcut; gz<gcut+1; ++gz){

        if(gx==0 && gy==0 && gz==0) continue;
            
        Eigen::Matrix<double,1,3> gg {1.0*gx, 1.0*gy, 1.0*gz};

        Eigen::Matrix<double,1,3> gvec = gg * gvecs;
        auto gdiff2 = gvec.dot(gvec);

        for (auto jj=0; jj<nq; ++jj){
          auto gr = gvec.dot(ts[ia]-ts[jj]);
          longrange = longrange + 4.0*M_PI/vol * gvec.transpose()/gdiff2 * sin(gr) * exp(-0.5*sig*sig * gdiff2) * qs[jj];
        }
      }
    }
  }


  return shortrange+longrange;
}



#endif
