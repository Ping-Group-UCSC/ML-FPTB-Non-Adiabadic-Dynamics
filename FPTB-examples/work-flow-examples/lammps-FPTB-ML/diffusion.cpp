#include "perovskite_tb.cpp"
#include <math.h>
#include <cmath>
#include "blot/blotmath.hpp"
#include "blot/units.hpp"
#include "auxparam.hpp"
#include <mpi.h>


void set_val(int ind, std::vector<double> hop, MatWrap& dx, MatWrap& dy, MatWrap& dz, int row_start_x, int row_end_x, int row_start_y, int row_end_y, int row_start_z, int row_end_z){


  if ( ind-1 > row_start_x-1 && ind-1 < row_end_x ) {
    MatSetValue(dx.A, ind-1, ind-1, hop[0] , ADD_VALUES);
  }

  if ( ind-1 > row_start_y-1 && ind-1 < row_end_y ) {
    MatSetValue(dy.A, ind-1, ind-1, hop[1], ADD_VALUES);
  }

  if ( ind-1 > row_start_z-1 && ind-1 < row_end_z ) {
    MatSetValue(dz.A, ind-1, ind-1, hop[2] , ADD_VALUES);
  }
}


void evaluateDr(VecWrap& psi, MatWrap& r, PetscScalar& Pdr2P_val){

  VecWrap rPsi(psi.n);
  PetscScalar PrP_val;

  MatMult(r.A, psi.V, rPsi.V);

  VecDot(psi.V, rPsi.V, &PrP_val);
  Mat PrPMat;
  MatCreateConstantDiagonal(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, psi.n, psi.n, PrP_val, &PrPMat);

  MatAXPY(r.A, -1.0, PrPMat,SAME_NONZERO_PATTERN); 

  Mat dr2;
  MatMatMult(r.A, r.A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dr2);
  VecWrap dr2Psi(psi.n);

  MatMult(dr2, psi.V, dr2Psi.V);


  VecDot(psi.V, dr2Psi.V, &Pdr2P_val);
}

//calculate diffusion
void get_diffusion(PerovData data, std::vector<Eigen::Vector3d>& coord, VecWrap& psi){

  MatWrap dx(data.norb,data.norb);
  MatWrap dy(data.norb,data.norb);
  MatWrap dz(data.norb,data.norb);
  
  int row_start_x, row_end_x;
  MatGetOwnershipRange(dx.A, &row_start_x, &row_end_x);
  int row_start_y, row_end_y;
  MatGetOwnershipRange(dy.A, &row_start_y, &row_end_y);
  int row_start_z, row_end_z;
  MatGetOwnershipRange(dz.A, &row_start_z, &row_end_z);

  for (int iatom=0; iatom<data.nn;++iatom){

    Atomehr atom1(iatom);
    double xc, yc, zc;
    vector<int> psiInd(1);
  
    xc = coord[iatom](0);
    yc = coord[iatom](1);
    zc = coord[iatom](2);
    for (auto orb1 : data.orblut.lookup_idx(Atomehr(atom1))){

        std::vector<double> pos = { xc, yc, zc};

        set_val(orb1.value, pos, dx, dy, dz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

    }

  }


  MatAssemblyBegin(dx.A, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(dy.A, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(dz.A, MAT_FINAL_ASSEMBLY);

  MatAssemblyEnd(dx.A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dy.A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dz.A, MAT_FINAL_ASSEMBLY);

  PetscScalar Dx,Dy,Dz;

  evaluateDr(psi, dx, Dx);
  evaluateDr(psi, dy, Dy);
  evaluateDr(psi, dz, Dz);

  

  if (me==0){ 
    std::ofstream out;
    out.open("diff.xyz",std::ios::app);
    out << std::to_string(abs(Dx)) << " " << std::to_string(abs(Dy)) << " " << std::to_string(abs(Dz)) << " " << endl;
    out.close();
  }
}

