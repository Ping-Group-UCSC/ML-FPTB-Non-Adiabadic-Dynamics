#ifndef PEROVSKITE_TB
#define PEROVSKITE_TB

#include "blot/lookuptable.hpp"
#include "blot/tight_binding.hpp"
#include "Eigen/Dense"
#include <algorithm>
#include <cmath>

class Atom_lut{};
class Atom_typ_lut{};
class Atom_ix_lut{};
class Atom_iy_lut{};
class Atom_iz_lut{};

class Orb_lut{};
class Orb_typ_lut{};
class Orb_ix_lut{};
class Orb_iy_lut{};
class Orb_iz_lut{};
class Ddir_lut{};
class Pdir_lut{};
class Spin_lut{};

class Shell_idx_lut{};
class Shell_lut{};
class Atom_nbr_lut{};
class Dx_lut{};
class Dy_lut{};
class Dz_lut{};

//class Hop_idx_lut{};
//class Hop_lut{};
//class Orb1_lut{};

typedef luttype<int,Atom_lut> Atom;
typedef luttype<string,Atom_typ_lut> Atom_typ;
typedef luttype<int,Atom_ix_lut> Atom_ix;
typedef luttype<int,Atom_iy_lut> Atom_iy;
typedef luttype<int,Atom_iz_lut> Atom_iz;
typedef lookuptable<Atom, Atom_typ, Atom_ix, Atom_iy, Atom_iz> Atomlut;

typedef luttype<int,Orb_lut> Orb;
typedef luttype<string,Orb_typ_lut> Orb_typ;
typedef luttype<int,Orb_ix_lut> Orb_ix;
typedef luttype<int,Orb_iy_lut> Orb_iy;
typedef luttype<int,Orb_iz_lut> Orb_iz;
typedef luttype<string,Ddir_lut> Ddir;
typedef luttype<string,Pdir_lut> Pdir;
typedef luttype<string,Spin_lut> Spin;
typedef lookuptable<Orb, Orb_typ, Atom, Ddir, Pdir, Spin> Orblut;

typedef luttype<int,Shell_idx_lut> Shell_idx;
typedef luttype<string,Shell_lut> Shell;
typedef luttype<int,Atom_nbr_lut> Atom_nbr;
typedef luttype<int,Dx_lut> Dx;
typedef luttype<int,Dy_lut> Dy;
typedef luttype<int,Dz_lut> Dz;
typedef lookuptable<Shell_idx, Shell, Atom, Atom_nbr, Dx, Dy, Dz> Shelllut;

//typedef luttype<int,Hop_idx_lut> Hop_idx;
//typedef luttype<string,Hop_lut> Hop;
//typedef luttype<int,Orb1_lut> Orb1;
//typedef lookuptable<Hop_idx,Hop,Orb,Orb1,Dx,Dy,Dz> Hoplut;

struct PerovData{
  Atomlut atomlut;
  Orblut orblut;
  Shelllut shelllut;
  int nx, ny, nz;//number of unit cells
  int nn; //number of atoms
  Eigen::MatrixXd cell;
  std::vector<double> qs; //charges

  PerovData(Atomlut atomlut1,
      Orblut orblut1,
      Shelllut shelllut1,
      int nx1, int ny1, int nz1, int nn1,
      const Eigen::MatrixXd & cell1);
};


PerovData::PerovData(Atomlut atomlut1,
    Orblut orblut1,
    Shelllut shelllut1,
    int nx1, int ny1, int nz1, int nn1,
    const Eigen::MatrixXd & cell1) :
  atomlut(atomlut1), orblut(orblut1), shelllut(shelllut1), nx(nx1), ny(ny1), nz(nz1), nn(nn1), cell(cell1) {
  
  //initialize charges  
    for(auto ia=0; ia<nn; ++ia){
      auto atom_typ = atomlut.lookup_col<Atom_typ>(Atom(ia));
      if (atom_typ == "Pb"){
        qs.push_back(2.0);
      }
      else if (atom_typ == "Br1" || atom_typ == "Br2" || atom_typ == "Br3"){
        qs.push_back(-1.0);
      }
      else if (atom_typ == "Cs"){
        qs.push_back(1.0);
      }
    }
  }

template<class A>
Eigen::Vector3d xyzvec1(luttype<string,A> ddir){

  //std::transform(ddir.value.begin(), ddir.value.end(), ddir.value.begin(), (int (*)(int))tolower);
  if (ddir=="x"){
    Eigen::Vector3d ans {1.0,0,0};
    return ans;
  }
  else if (ddir=="y"){
    Eigen::Vector3d ans {0,1.0,0};
    return ans;
  }
  else if (ddir=="z"){
    Eigen::Vector3d ans {0,0,1.0};
    return ans;
  }

  Eigen::Vector3d z {0,0,0.0};
  return z;
}

template<class A>
Eigen::Vector3d xyzvec2(luttype<string,A> ddir){

  //std::transform(ddir.value.begin(), ddir.value.end(), ddir.value.begin(), (int (*)(int))tolower);
  if (ddir=="z"){
    Eigen::Vector3d ans {1.0,0,0};
    return ans;
  }
  else if (ddir=="x"){
    Eigen::Vector3d ans {0,1.0,0};
    return ans;
  }
  else if (ddir=="y"){
    Eigen::Vector3d ans {0,0,1.0};
    return ans;
  }

  Eigen::Vector3d z {0,0,0.0};
  return z;
}

Eigen::Vector3d posvec0(const Atomlut & atomlut, Atom atom1){
  //returns standard coordinates of an atom as a vector
  //  ix iy iz for Cs
  //  ix iy+0.5 iz+0.5 for Br1, etc
  //  ix+0.5 iy+0.5 iz+0.5 for Pb
  auto ix = atomlut.lookup_col<Atom_ix>(atom1);
  auto iy = atomlut.lookup_col<Atom_iy>(atom1);
  auto iz = atomlut.lookup_col<Atom_iz>(atom1);
  auto typ = atomlut.lookup_col<Atom_typ>(atom1);

  Eigen::Vector3d t0 {1.0*ix.value, 1.0*iy.value, 1.0*iz.value};
  
  if (typ=="Pb") {
    Eigen::Vector3d t1 {0.5,0.5,0.5};
    return t0+t1;
  }
  else if (typ=="Br1") {
    Eigen::Vector3d t1 {0.0,0.5,0.5};
    return t0+t1;
  }
  else if (typ=="Br2") {
    Eigen::Vector3d t1 {0.5,0.0,0.5};
    return t0+t1;
  }
  else if (typ=="Br3") {
    Eigen::Vector3d t1 {0.5,0.5,0.0};
    return t0+t1;
  }
  else if (typ=="Cs") {
    Eigen::Vector3d t1 {0.0,0.0,0.0};
    return t0+t1;
  }

  return t0;
}

//translate atoms to be centered on standard positions:
// LX0,LY0,LZ0 = unitcell dimensions
Eigen::Vector3d translate_posvec0(const Atomlut & atomlut, Atom atom, const Eigen::Vector3d & c_old, double LX0, double LY0, double LZ0){
  Eigen::Vector3d pv = posvec0(atomlut,atom) ;
  Eigen::Vector3d uc {LX0, LY0, LZ0};
  Eigen::Vector3d uc1 {pv(0)*LX0, pv(1)*LY0, pv(2)*LZ0};
  Eigen::Vector3d c_new = c_old - uc1 + 0.5*uc;
  for (auto ii=0; ii<3; ++ii){
    c_new(ii) = fmod(c_new(ii), uc(ii));
  }
  c_new = c_new + uc1 - 0.5*uc;
  return c_new;
}

//
// FrameFn
//

struct pbp_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));

    auto axis1 = xyzvec1(pdir1);
    auto axis2 = xyzvec2(pdir1);
    auto axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct pbs_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    return Eigen::MatrixXd::Identity(3,3);
  }
};

struct brp_par_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));

    auto axis1 = xyzvec1(pdir1);
    auto axis2 = xyzvec2(pdir1);
    auto axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct brp_perp_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto ddir1 = data.orblut.lookup_col<Ddir>(Orb(orb1));
    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));

    Eigen::Vector3d axis1 = xyzvec1(pdir1);
    Eigen::Vector3d axis2 = xyzvec2(ddir1);
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct pppi_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));

    Eigen::Vector3d axis1 = xyzvec1(pdir1);
    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d axis2 = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    axis2.normalize();
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct ppsig_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir = data.orblut.lookup_col<Pdir>(Orb(orb1));

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));

    Eigen::Vector3d axis1 = xyzvec2(pdir);
    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d axis2 = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    axis2.normalize();
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct spsig_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));
    auto atom_typ1 = data.atomlut.lookup_col<Atom_typ>(atom1);

    Pdir pdir;
    if (atom_typ1 == "Pb"){
      pdir = data.orblut.lookup_col<Pdir>(Orb(orb2));
    }
    else{
      pdir = data.orblut.lookup_col<Pdir>(Orb(orb1));
    }

    Eigen::Vector3d axis1 = xyzvec2(pdir);
    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d axis2 = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    axis2.normalize();
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

bool abs_comparator(double a, double b){
  return abs(a) < abs(b);
}

struct ppn_br_to_pb_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));


    Eigen::Vector3d axis1 = xyzvec1(pdir2);
    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d axis2 = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    axis2.normalize();

    auto v = std::max_element(axis2.begin(), axis2.end(), abs_comparator);
    if (signbit(*v)) axis1 = axis1 * (-1.0); 

    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct ppn_pb_to_br_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));


    Eigen::Vector3d axis1 = xyzvec1(pdir2);
    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d axis2 = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    axis2.normalize();

    auto v = std::max_element(axis2.begin(), axis2.end(), abs_comparator);
    if (signbit(*v)) axis1 = axis1 * (-1.0); 

    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct spn_pb_to_br_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));
    auto atom_typ1 = data.atomlut.lookup_col<Atom_typ>(atom1);

    Pdir pdir;
    if (atom_typ1 == "Pb"){
      pdir = data.orblut.lookup_col<Pdir>(Orb(orb2));
    }
    else{
      pdir = data.orblut.lookup_col<Pdir>(Orb(orb1));
    }

    Eigen::Vector3d axis1 = xyzvec1(pdir);
    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d axis2 = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    axis2.normalize();
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct socpb_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));

    Eigen::Vector3d axis1 = xyzvec1(pdir1);
    Eigen::Vector3d axis2 = xyzvec1(pdir2);
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};



struct socbr_par_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));
    auto ddir1 = data.orblut.lookup_col<Ddir>(Orb(orb1));
    auto ddir2 = data.orblut.lookup_col<Ddir>(Orb(orb2));

    Eigen::Vector3d axis1;
    Eigen::Vector3d axis2;
    if (ddir1 == pdir1){
      axis1 = xyzvec1(pdir1);
      axis2 = xyzvec1(pdir2);
    }
    else if (ddir2 == pdir2){
      axis1 = xyzvec1(pdir2);
      axis2 = xyzvec1(pdir1);
    }

    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

struct socbr_perp_FrameFn{
  Eigen::MatrixXd operator() (const PerovData & data, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));

    Eigen::Vector3d axis1 = xyzvec1(pdir1);
    Eigen::Vector3d axis2 = xyzvec1(pdir2);
    Eigen::Vector3d axis3 = axis1.cross(axis2);

    Eigen::MatrixXd ans(3,3);
    ans << axis1 , axis2, axis3;
    return ans;
  }
};

//
// FeatureFn
//

struct vcoul_FeatureFn{
  std::vector<double> operator() (const PerovData & data, const Eigen::MatrixXd & frame, std::vector<int> hop,  const std::vector<Eigen::Vector3d> & atompos){
    int orb1 = (int)hop[3];
    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto vcoul = ewald_phi(atom1.value, data.cell, atompos, data.qs, 12.0,2,2);
    std::vector<double> ans {vcoul};


    return ans;
  }
};

struct xyz_FeatureFn{
  std::vector<double> operator() (const PerovData & data, const Eigen::MatrixXd & frame, std::vector<int> hop, const std::vector<Eigen::Vector3d> & atompos){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));

    Eigen::Vector3d dd {1.0*dx, 1.0*dy, 1.0*dz};
    Eigen::Vector3d rdiff = atompos[atom2.value] - atompos[atom1.value] + data.cell*dd;
    auto xdiff = rdiff.dot(frame(Eigen::all,1));//axis 2 was defined in the the frames to be the bonding direction
    auto ydiff = rdiff.dot(frame(Eigen::all,0));
    auto zdiff = rdiff.dot(frame(Eigen::all,2));

    
    std::vector<double> ans {xdiff, ydiff, zdiff};

    return ans;
  }
};

//
// ModelFn
//

struct vcoul_ModelFn{
  std::complex<double> operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double vcoul = feats[0];
    return params[0]*vcoul + params[1];
  }
};

struct spsig_ModelFn{
  std::complex<double> operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    return params[0]*exp(params[1]*x + params[2]*(y*y + z*z)) + params[3];
  }
};

struct ppsig_ModelFn{
  std::complex<double> operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    return params[0]*exp(params[1]*x + params[2]*(y*y + z*z)) + params[3];
  }
};

struct pppi_ModelFn{
  std::complex<double> operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    return params[0]*exp(params[1]*x + params[2]*y*y + params[3]*z*z) + params[4];
  }
};

struct nonbonding_ModelFn{
  std::complex<double> operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    return params[0]*y*exp(params[1]*x + params[2]*z*z);
  }
};

//
// PhaseFn
//

struct spsig_PhaseFn{
  std::complex<double> operator() (const PerovData & data, std::vector<int> hop){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atom>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atom>(Orb(orb2));

    Eigen::Vector3d nd {1.0*data.nx*dx, 1.0*data.ny*dy, 1.0*data.nz*dz};
    Eigen::Vector3d dr = posvec0(data.atomlut,atom2) - posvec0(data.atomlut,atom1) + nd;
    if (dr(0)>0.25 || dr(1)>0.25 || dr(2)>0.25 ){
      return 1.0;
    }
    else{
      return -1.0;
    }
  }
};

struct id_PhaseFn{
  std::complex<double> operator() (const PerovData & data, std::vector<int> hop){
    return 1.0;
  }
};

struct soc_PhaseFn{
  std::complex<double> operator() (const PerovData & data, std::vector<int> hop){
    int dx = (int)hop[0];
    int dy = (int)hop[1];
    int dz = (int)hop[2];
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto spin1 = data.orblut.lookup_col<Spin>(Orb(orb1));
    auto spin2 = data.orblut.lookup_col<Spin>(Orb(orb2));
    auto pdir1 = data.orblut.lookup_col<Pdir>(Orb(orb1));
    auto pdir2 = data.orblut.lookup_col<Pdir>(Orb(orb2));

    if(spin1=="up" && pdir1=="x" && spin2=="up" && pdir2=="y") return 1.0i; 
    if(spin1=="dn" && pdir1=="x" && spin2=="dn" && pdir2=="y") return -1.0i; 
    if(spin1=="dn" && pdir1=="x" && spin2=="up" && pdir2=="z") return -1.0; 
    if(spin1=="up" && pdir1=="x" && spin2=="dn" && pdir2=="z") return 1.0; 
    if(spin1=="up" && pdir1=="y" && spin2=="dn" && pdir2=="z") return 1.0i; 
    if(spin1=="dn" && pdir1=="y" && spin2=="up" && pdir2=="z") return 1.0i; 

    return 0.0;

  }
};

template <class HM>
void insert_hop_amplitudes(std::vector<std::vector<double>> & hop_amplitudes, HM hm, std::vector<std::vector<int>> hops, std::vector<double> params, std::vector<Eigen::Vector3d> coord){
  for (auto hop : hops){

    complex<double> hval = hm(hop, params, coord);
    std::vector<double> hop_amplitude {1.0*hop[0], 1.0*hop[1], 1.0*hop[2], 1.0*hop[3], 1.0*hop[4], hval.real(), hval.imag()};
    hop_amplitudes.push_back(hop_amplitude);
  }
};



#endif
