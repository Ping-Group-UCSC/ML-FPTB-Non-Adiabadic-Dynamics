#ifndef AUXPARAM
#define AUXPARAM

#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "domain.h"
#include "library.h"
#include <Eigen/Dense>


using namespace LAMMPS_NS;

  int me,nprocs;
  int is,ie;
  int iscd,iecd;

  LAMMPS *lmp;

  int numhop,icounthop,jcounthop;

  int *ndatahopamp;
  int *ndisplhoparam;

  int *ndatacd;
  int *displcd;

  Eigen::Matrix3d gvecs;
  Eigen::Matrix3d cell;

  double vol;

  double *tewald_phi;

  const double sigewald=12.0;

  double sqrt2M_PIsig,fM_PIvol,hsigsig,onesqrt2sig,hpsigpsig,coefewald_field;

#endif
