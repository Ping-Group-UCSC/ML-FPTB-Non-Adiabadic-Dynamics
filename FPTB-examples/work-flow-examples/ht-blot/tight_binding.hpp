#ifndef TIGHT_BINDING
#define TIGHT_BINDING

#include <slepceps.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include "blot/blotio.hpp"
#include "blot/blotmath.hpp"
#include "blot/basics.hpp"


//class for making a tight binding hamiltonian from hoppings
struct TBHammaker{
  int dim; //number of spatial dimensions
  int nhs; // size of hilbert space
  std::vector<std::vector<double> > & hopsdat;

  PetscErrorCode create(MatWrap & ans, std::vector<double> const & kpt){

    PetscErrorCode ierr;

    for( auto row : hopsdat){
      double kd = 0.0; // kpt dot dcell
      for(auto ii=0; ii!=dim; ii++){
          kd = kd + row[ii]*kpt[ii];
      }

      int from = row[dim] - 1; //the hopsdat is 1-indexed
      int to = row[dim+1] - 1;
      PetscScalar val = row[dim+2] + row[dim+3]*PETSC_i;

      auto phase = PetscExpComplex(2.0 * PETSC_PI * PETSC_i * kd);

      ierr = MatSetValue(ans.A,from,to,phase*val,ADD_VALUES);CHKERRQ(ierr); 

      if (any_of(row.begin(), row.begin()+dim, [](double j){return abs(j)>0.5;} ) || from!=to) {
        ierr = MatSetValue(ans.A,to,from,PetscConjComplex(phase*val),ADD_VALUES);CHKERRQ(ierr); 
      }

    }

    ierr = MatAssemblyBegin(ans.A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(ans.A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    return ierr;
  }
};

struct Solve_ham_default_options{
  PetscErrorCode apply_options(EPS eps){
    PetscErrorCode ierr=0;
    return ierr;
  }
};

struct Solve_ham_dense_options{
  PetscErrorCode apply_options(EPS eps){
    PetscErrorCode ierr=0;
    ierr = EPSSetType(eps, EPSLAPACK);
    return ierr;
  }
};

struct Solve_ham_interval_options{
  PetscReal a;
  PetscReal b;
  PetscErrorCode apply_options(EPS eps){
    PetscErrorCode ierr=0;
    ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);
    ierr = EPSSetInterval(eps,a,b);
    return ierr;
  }
};

struct Solve_ham_target_options{
  PetscReal a;
  PetscInt nev;
  PetscErrorCode apply_options(EPS eps){
    PetscErrorCode ierr=0;
    ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_REAL);
    ierr = EPSSetTarget(eps,a);
    //ierr = EPSSetDimensions(eps, nev, 2*nev, nev);
    ierr = EPSSetDimensions(eps, nev, PETSC_DEFAULT, PETSC_DEFAULT);
    return ierr;
  }
};



//Driver for solving TB hamiltonian
//Htype contains the information for constructing the Hamiltonian:
//  PetscErrorCode Htype::create(MatWrap & ham, Ptype param);
//  int Htype::nhs
//Ptype are the parameters  
//Options specify how to carry out the diagonalization
//  PetscErrorCode Options::apply_options(EPS eps);
template <class Htype, class Ptype, class Options>
PetscErrorCode solve_ham_params(Htype hamfun, std::vector<Ptype> params, bool save_wfns, Options op){

  PetscErrorCode ierr;

  //ierr = SlepcInitWrap();

  for (int iparam=0; iparam!=params.size(); ++iparam){
    cout << "doing iparam = " << iparam << endl;
    cout << flush;
    MatWrap ham(hamfun.nhs,hamfun.nhs);
    EPS eps;

    auto param = params[iparam];
    ierr = hamfun.create(ham, param);

    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
    ierr = EPSSetOperators(eps,ham.A,NULL);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

    ierr = op.apply_options(eps);CHKERRQ(ierr);
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

    ierr = EPSSolve(eps);CHKERRQ(ierr);

    PetscInt nconv;
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);

    mkdir("work",0777);

    std::vector<PetscReal> evals; 
    PetscReal error;
    PetscScalar evalr,evali;

    if (nconv>0) {

      //get eigenvalues first
      for (auto i=0;i!=nconv;i++) {
        ierr = EPSGetEigenpair(eps,i,&evalr,&evali,NULL,NULL);CHKERRQ(ierr);
        ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
        evals.push_back(PetscRealPart(evalr));
      }

      auto sortidxs = sortperm(evals);
      sort(evals.begin(), evals.end());

      //get eigenvectors if requested
      if(save_wfns){
        mkdir(("work/"+to_string(iparam)).c_str(), 0777);
        for (auto ieig=0; ieig!=nconv; ieig++) {

          auto i = sortidxs[ieig];

          Vec xr, xi;
          ierr = EPSGetEigenpair(eps,i,&evalr,&evali,xr,xi);CHKERRQ(ierr);

          //storing eigenvectors
          savetxt("work/"+to_string(iparam)+"/evec"+to_string(ieig)+"_real.dat", xr);
          savetxt("work/"+to_string(iparam)+"/evec"+to_string(ieig)+"_imag.dat", xi);

          ierr = VecDestroy(&xr);CHKERRQ(ierr);
          ierr = VecDestroy(&xi);CHKERRQ(ierr);
        }
      }
    }
      
    mkdir(("work/"+to_string(iparam)).c_str(), 0777);
    savetxt("work/"+to_string(iparam)+"/evals.dat",evals);


    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  }

  ierr = SlepcFinalize(); 

  return ierr;
}

//template for a functor that computes hopping amplitudes given atomic positions
//
//SystemData: contains properties like unitcell size, number of atoms, etc
//FrameFn: a functor that assigns a coordinate frame for each hopping
//FeatureFn: a functor that computes features from atomic positions
//ModelFn: a functor that computes hopping amplitudes from features
//PhaseFn: a functor that performs a phase transformation on the hopping amplitudes
//
//vector<int> hop {dx, dy, dz, orb1, orb2};
template<class SystemData, class FrameFn, class FeatureFn, class ModelFn, class PhaseFn >
struct HoppingModel{
  SystemData data;
  FrameFn framefn;
  FeatureFn featurefn;
  ModelFn modelfn;
  PhaseFn phasefn;

  HoppingModel(SystemData data1) : data(data1) { };

  complex<double> operator() (std::vector<int> hop, std::vector<double> params, std::vector<Eigen::Vector3d> atompos){
    auto frame = framefn(data,hop,atompos);
    auto feats = featurefn(data,frame,hop,atompos);
    auto hval = modelfn(data,params,feats);
    auto hval1 = hval / phasefn(data,hop);

    return hval1;
  }
};



#endif
