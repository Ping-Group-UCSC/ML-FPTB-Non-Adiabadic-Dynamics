#include "perovskite_tb.cpp"
#include <math.h>
#include <cmath>
#include "blot/blotmath.hpp"
#include "blot/units.hpp"
#include "auxparam.hpp"
#include <mpi.h>

//calculate charge density 
// [coord] = angstrom
// [md_charges] = e
void save_charge(PerovData data, std::vector<Eigen::Vector3d>& coord, std::vector<double> md_charges, VecWrap& psi){

  std::vector<double> charge_density(data.nn, 0.0);

  for (auto atom : data.atomlut.lookup_idx()){
    std::vector<int> iorbs ;
    for (auto iorb : data.orblut.lookup_idx(Atomehr(atom))){
        iorbs.push_back(iorb.value-1);
    }

    VecWrap temp = psi.get(iorbs);

    PetscScalar cdens;
    VecDot(temp.V, temp.V, &cdens);

    charge_density[atom.value] = charge_density[atom.value] + abs(cdens);

  }
  if (me==0){
  
    double chg_tot=0.0;
    
    std::ofstream out;
    out.open("charge.xyz",std::ios::app);
    std::ofstream outC;
    outC.open("total_charge.dat",std::ios::app);
    out << "step" << endl;
    for (int ii=0; ii<data.nn; ++ii){
      chg_tot += charge_density[ii];
      out << charge_density[ii] << endl;
      }
    outC << chg_tot << endl;
    out.close();
    outC.close();
  }

}

