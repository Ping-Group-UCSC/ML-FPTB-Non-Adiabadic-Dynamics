#ifndef PEROVSKITE_TB
#define PEROVSKITE_TB

#include "blot/lookuptable.hpp"
#include "blot/tight_binding.hpp"
#include "Eigen/Dense"
#include <algorithm>
#include <cmath>

#include "auxparam.hpp"

class Atomehr_lut{};
class Atomehr_typ_lut{};
class Atomehr_ix_lut{};
class Atomehr_iy_lut{};
class Atomehr_iz_lut{};

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
class Atomehr_nbr_lut{};
class Dx_lut{};
class Dy_lut{};
class Dz_lut{};

//class Hop_idx_lut{};
//class Hop_lut{};
//class Orb1_lut{};

typedef luttype<int,Atomehr_lut> Atomehr;
typedef luttype<string,Atomehr_typ_lut> Atomehr_typ;
typedef luttype<int,Atomehr_ix_lut> Atomehr_ix;
typedef luttype<int,Atomehr_iy_lut> Atomehr_iy;
typedef luttype<int,Atomehr_iz_lut> Atomehr_iz;
typedef lookuptable<Atomehr, Atomehr_typ, Atomehr_ix, Atomehr_iy, Atomehr_iz> Atomehrlut;

typedef luttype<int,Orb_lut> Orb;
typedef luttype<string,Orb_typ_lut> Orb_typ;
typedef luttype<int,Orb_ix_lut> Orb_ix;
typedef luttype<int,Orb_iy_lut> Orb_iy;
typedef luttype<int,Orb_iz_lut> Orb_iz;
typedef luttype<string,Ddir_lut> Ddir;
typedef luttype<string,Pdir_lut> Pdir;
typedef luttype<string,Spin_lut> Spin;
typedef lookuptable<Orb, Orb_typ, Atomehr, Ddir, Pdir, Spin> Orblut;

typedef luttype<int,Shell_idx_lut> Shell_idx;
typedef luttype<string,Shell_lut> Shell;
typedef luttype<int,Atomehr_nbr_lut> Atomehr_nbr;
typedef luttype<int,Dx_lut> Dx;
typedef luttype<int,Dy_lut> Dy;
typedef luttype<int,Dz_lut> Dz;
typedef lookuptable<Shell_idx, Shell, Atomehr, Atomehr_nbr, Dx, Dy, Dz> Shelllut;

//typedef luttype<int,Hop_idx_lut> Hop_idx;
//typedef luttype<string,Hop_lut> Hop;
//typedef luttype<int,Orb1_lut> Orb1;
//typedef lookuptable<Hop_idx,Hop,Orb,Orb1,Dx,Dy,Dz> Hoplut;

typedef std::vector<std::vector<std::vector<int>>> Vcutlut;

struct PerovData{
  Atomehrlut atomlut;
  Orblut orblut;
  Shelllut shelllut;
  Vcutlut vcutlut;
  int nx, ny, nz;//number of unit cells
  int nn; //number of atoms
  int norb; //number of orbitals
  Eigen::MatrixXd cell;
  std::vector<double> qs; //charges

  std::vector<double> pbp_params;
  std::vector<double> pbs_params ;
  std::vector<double>  brp_par_params ;
  std::vector<double> brp_perp_params ;
  std::vector<double>  spsig_params ;
  std::vector<double>  ppsig_params ;
  std::vector<double>   pppi_params ;
  std::vector<double>  ppn_br_to_pb_params ;
  std::vector<double> ppn_pb_to_br_params ;
  std::vector<double>  spn_pb_to_br_params ;
  std::vector<double> socpb_params ;
  std::vector<double> socbr_par_params ;
  std::vector<double>  socbr_perp_params ;


  PerovData(Atomehrlut atomlut1,
      Orblut orblut1,
      Shelllut shelllut1,
      Vcutlut vcutlut1,
      int nx1, int ny1, int nz1, int nn1, int norb1,
      const Eigen::MatrixXd & cell1
      );

  int pbp_idx( Atomehr ia, Spin sp, Pdir pdir) const;
  int brp_idx( Atomehr ia, Spin sp, Pdir pdir) const;
};


PerovData::PerovData(Atomehrlut atomlut1,
    Orblut orblut1,
    Shelllut shelllut1,
    Vcutlut vcutlut1,
    int nx1, int ny1, int nz1, int nn1, int norb1,
    const Eigen::MatrixXd & cell1
    ) :
  atomlut(atomlut1), orblut(orblut1), shelllut(shelllut1), vcutlut(vcutlut1), nx(nx1), ny(ny1), nz(nz1), nn(nn1), norb(norb1), cell(cell1) {
  
  //initialize charges  
    for(auto ia=0; ia<nn; ++ia){
      auto atom_typ = atomlut.lookup_col<Atomehr_typ>(Atomehr(ia));
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

    //load parameters
    pbp_params = loadtxtvec<double>("data_to_load/pbp_params.dat");
    pbs_params = loadtxtvec<double>("data_to_load/pbs_params.dat");
    brp_par_params = loadtxtvec<double>("data_to_load/brp_par_params.dat");
    brp_perp_params = loadtxtvec<double>("data_to_load/brp_perp_params.dat");
    spsig_params = loadtxtvec<double>("data_to_load/spsig_params.dat");
    ppsig_params = loadtxtvec<double>("data_to_load/ppsig_params.dat");
    pppi_params = loadtxtvec<double>("data_to_load/pppi_params.dat");
    ppn_br_to_pb_params = loadtxtvec<double>("data_to_load/ppn_br_to_pb_params.dat");
    ppn_pb_to_br_params = loadtxtvec<double>("data_to_load/ppn_pb_to_br_params.dat");
    spn_pb_to_br_params = loadtxtvec<double>("data_to_load/spn_pb_to_br_params.dat");
    socpb_params = loadtxtvec<double>("data_to_load/socpb_params.dat");
    socbr_par_params = loadtxtvec<double>("data_to_load/socbr_par_params.dat");
    socbr_perp_params = loadtxtvec<double>("data_to_load/socbr_perp_params.dat");

  }


//WARNING: orb index hard-coded here must be consistent with orb_index.dat!
int PerovData::pbp_idx( Atomehr ia, Spin sp, Pdir pdir) const{

  int npbs = 2*nx*ny*nz;
  int npb = nx*ny*nz;
  int ncs = nx*ny*nz;
  int pdir_idx;
  if (pdir=="z"){
    pdir_idx = 0;
  }
  else if (pdir=="x"){
    pdir_idx = 1;
  }
  else if (pdir=="y"){
    pdir_idx = 2;
  }

  int sp_idx;
  if (sp=="up"){
    sp_idx = 0;
  }
  else if (sp=="dn"){
    sp_idx = 1;
  }

  return npbs + 6*(ia.value-ncs) + 2*pdir_idx + sp_idx +1;
}

//WARNING: orb index hard-coded here must be consistent with orb_index.dat!
int PerovData::brp_idx( Atomehr ia, Spin sp, Pdir pdir) const{

  int npbs = 2*nx*ny*nz;
  int npbp = 6*nx*ny*nz;
  int npb = nx*ny*nz;
  int ncs = nx*ny*nz;
  int pdir_idx;
  if (pdir=="z"){
    pdir_idx = 0;
  }
  else if (pdir=="x"){
    pdir_idx = 1;
  }
  else if (pdir=="y"){
    pdir_idx = 2;
  }

  int sp_idx;
  if (sp=="up"){
    sp_idx = 0;
  }
  else if (sp=="dn"){
    sp_idx = 1;
  }

  return npbs +npbp + 6*(ia.value - ncs- npb) + 2*pdir_idx + sp_idx +1;
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

Eigen::Vector3d posvec0(const Atomehrlut & atomlut, Atomehr atom1){
  //returns standard coordinates of an atom as a vector  //  ix iy iz for Cs
  //  ix iy iz for Cs
  //  ix iy+0.5 iz+0.5 for Br1, etc
  //  ix+0.5 iy+0.5 iz+0.5 for Pb
  auto ix = atomlut.lookup_col<Atomehr_ix>(atom1);
  auto iy = atomlut.lookup_col<Atomehr_iy>(atom1);
  auto iz = atomlut.lookup_col<Atomehr_iz>(atom1);
  auto typ = atomlut.lookup_col<Atomehr_typ>(atom1);

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
Eigen::Vector3d translate_posvec0(const Atomehrlut & atomlut, Atomehr atom, const Eigen::Vector3d & c_old, double LX0, double LY0, double LZ0){
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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));
    auto atom_typ1 = data.atomlut.lookup_col<Atomehr_typ>(atom1);

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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));
    auto atom_typ1 = data.atomlut.lookup_col<Atomehr_typ>(atom1);

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
// AG Cutoff Vcoul Function
//


double coul_cutoff_AG(const PerovData & data, int a1, const std::vector<Eigen::Vector3d> & atompos, int vcut_x, int vcut_y, int vcut_z){

  //int orb1 = (int)hop[3];
  //auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
  //auto atom1 = data.atomlut.lookup_col<Atomehr>(a1) ;
  auto atom1 = Atomehr(a1) ;
  auto ix = data.atomlut.lookup_col<Atomehr_ix>(atom1) ;
  auto iy = data.atomlut.lookup_col<Atomehr_iy>(atom1) ;
  auto iz = data.atomlut.lookup_col<Atomehr_iz>(atom1) ;

  double vcoul = 0.0;
  for (auto nbr : data.vcutlut[atom1.value]){
    auto dx = nbr[2];
    auto dy = nbr[3];
    auto dz = nbr[4];
    auto atom2 = nbr[1];

    Eigen::Matrix<double,3,1> dd {1.0*dx, 1.0*dy, 1.0*dz};

    double q2 = data.qs[atom2];
    auto rr = atompos[atom2] - atompos[atom1.value] + data.cell*dd;
    auto rdiff = rr.norm();
    vcoul = vcoul + q2/rdiff;

    }

  //std::vector<double> ans {vcoul};
  //return ans;
  return vcoul;
}



//
// FeatureFn
//

struct vcoul_FeatureFn{
  std::vector<double> operator() (const PerovData & data, const Eigen::MatrixXd & frame, std::vector<int> hop,  const std::vector<Eigen::Vector3d> & atompos){
    int orb1 = (int)hop[3];
    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    //auto vcoul = ewald_phi(atom1.value, data.cell, atompos, data.qs, 12.0,2,2);
    double vcoul=tewald_phi[atom1.value];
    std::vector<double> ans {vcoul};


    return ans;
  }
};

struct vcoul_cutoff_FeatureFn{
  std::vector<double> operator() (const PerovData & data, const Eigen::MatrixXd & frame, std::vector<int> hop,  const std::vector<Eigen::Vector3d> & atompos){
    int orb1 = (int)hop[3];
    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto ix = data.atomlut.lookup_col<Atomehr_ix>(atom1) ;
    auto iy = data.atomlut.lookup_col<Atomehr_iy>(atom1) ;
    auto iz = data.atomlut.lookup_col<Atomehr_iz>(atom1) ;

    double vcoul = 0.0;
    for (auto nbr : data.vcutlut[atom1.value]){
      auto dx = nbr[2];
      auto dy = nbr[3];
      auto dz = nbr[4];
      auto atom2 = nbr[1];

      Eigen::Matrix<double,3,1> dd {1.0*dx, 1.0*dy, 1.0*dz};

      double q2 = data.qs[atom2];
      auto rr = atompos[atom2] - atompos[atom1.value] + data.cell*dd;
      auto rdiff = rr.norm();
      vcoul = vcoul + q2/rdiff;

    }

    std::vector<double> ans {vcoul};
    return ans;
  }
};
//derivative of the coulomb potential at ia, with respect to the postion of atom iad
struct fcoul_cutoff_FeatureFn{
  std::vector<Eigen::RowVector3d> operator() (const PerovData & data, const std::vector<Eigen::Vector3d> & atompos, Atomehr ia, Atomehr iad, Dx dx, Dy dy, Dz dz){
    double q1 = data.qs[iad.value];
    Eigen::Matrix<double,3,1> dd {1.0*dx.value, 1.0*dy.value, 1.0*dz.value};
    auto rr = atompos[iad.value] - atompos[ia.value] - data.cell*dd;
    auto rdiff = rr.norm();

    double xderiv = -1.0*q1/pow(rdiff,3) * rr(0);
    double yderiv = -1.0*q1/pow(rdiff,3) * rr(1);
    double zderiv = -1.0*q1/pow(rdiff,3) * rr(2);

    Eigen::RowVector3d xyzderiv;
    xyzderiv << xderiv , yderiv, zderiv;
    std::vector<Eigen::RowVector3d> ans {xyzderiv};
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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

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
// DFeatureFn
//

//this are the gradient of the features
//ans(i,j) = d(feature_j)/d(direction_i)


//derivative of the coulomb potential at hop, with respect to the postion of atom iad
struct vcoul_cutoff_DFeatureFn{
  Eigen::MatrixXd  operator() (const PerovData & data, const std::vector<Eigen::Vector3d> & atompos, Atomehr ia, Atomehr iad, Dx dx, Dy dy, Dz dz){
  //Eigen::MatrixXd  operator() (const PerovData & data, const Eigen::MatrixXd & frame, std::vector<int> hop,  const std::vector<Eigen::Vector3d> & atompos, Atom iad, Dx dx, Dy dy, Dz dz){
    //int orb1 = (int)hop[3];
    //auto ia = data.orblut.lookup_col<Atom>(Orb(orb1));

    double xderiv = 0.0;
    double yderiv = 0.0;
    double zderiv = 0.0;
    if (ia == iad && dx==0 && dy==0 && dz==0){

      for (auto nbr : data.vcutlut[ia.value]){
        auto dx1 = nbr[2];
        auto dy1 = nbr[3];
        auto dz1 = nbr[4];
        auto atom2 = nbr[1];

        Eigen::Matrix<double,3,1> dd {1.0*dx1, 1.0*dy1, 1.0*dz1};

        double q2 = data.qs[atom2];
        auto rr = atompos[atom2] - atompos[ia.value] + data.cell*dd;
        auto rdiff = rr.norm();

        xderiv = xderiv +1.0*q2/pow(rdiff,3) * rr(0);
        yderiv = yderiv +1.0*q2/pow(rdiff,3) * rr(1);
        zderiv = zderiv +1.0*q2/pow(rdiff,3) * rr(2);
      }
    }

    else{ 
      double q1 = data.qs[iad.value];
      Eigen::Matrix<double,3,1> dd {1.0*dx.value, 1.0*dy.value, 1.0*dz.value};
      auto rr = atompos[iad.value] - atompos[ia.value] - data.cell*dd;
      auto rdiff = rr.norm();

      xderiv = -1.0*q1/pow(rdiff,3) * rr(0);
      yderiv = -1.0*q1/pow(rdiff,3) * rr(1);
      zderiv = -1.0*q1/pow(rdiff,3) * rr(2);
    }

    Eigen::Vector3d ans;
    ans << xderiv , yderiv, zderiv;
    return ans;
  }
};


struct xyz_DFeatureFn{
   Eigen::MatrixXd  operator() (const PerovData & data, const Eigen::MatrixXd & frame, std::vector<int> hop,  const std::vector<Eigen::Vector3d> & atompos, Atomehr iad, Dx dx, Dy dy, Dz dz){
    int orb1 = (int)hop[3];
    int orb2 = (int)hop[4];

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

    double sign = 0.0;

    if (atom1 == iad){ 
      sign = 1.0;
    }
    else if (atom2 == iad){
      sign = -1.0;
    }

    Eigen::RowVector3d perm1;
    perm1 << 0.0, 1.0, 0;
    Eigen::RowVector3d perm2;
    perm2 << 1.0, 0.0, 0;
    Eigen::RowVector3d perm3;
    perm3 << 0.0, 0, 1.0;

    Eigen::MatrixXd perm(3,3);
    perm << perm1 , perm2, perm3;

    return frame*perm*sign;
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

struct vcoul_deriv_ModelFn{
  Eigen::VectorXd operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    Eigen::VectorXd ans(1);
    ans << params[0];
    return ans;
  }
};

struct spsig_deriv_ModelFn{
  Eigen::Vector3d operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    double xderiv = params[0]*params[1]*exp(params[1]*x + params[2]*(y*y + z*z)) ;
    double yderiv = params[0]*2*y*params[2]*exp(params[1]*x + params[2]*(y*y + z*z)) ;
    double zderiv = params[0]*2*z*params[2]*exp(params[1]*x + params[2]*(y*y + z*z)) ;

    Eigen::Vector3d ans;
    ans << xderiv, yderiv, zderiv;
    return ans;
  }
};

struct ppsig_deriv_ModelFn{
  Eigen::Vector3d  operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    double xderiv = params[0]*params[1]*exp(params[1]*x + params[2]*(y*y + z*z)) ;
    double yderiv = params[0]*2*y*params[2]*exp(params[1]*x + params[2]*(y*y + z*z)) ;
    double zderiv = params[0]*2*z*params[2]*exp(params[1]*x + params[2]*(y*y + z*z)) ;

    Eigen::Vector3d ans;
    ans << xderiv, yderiv, zderiv;
    return ans;
  }
};

struct pppi_deriv_ModelFn{
  Eigen::Vector3d  operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    double xderiv = params[0]*params[1]*exp(params[1]*x + params[2]*y*y + params[3]*z*z) ;
    double yderiv = params[0]*2*y*params[2]*exp(params[1]*x + params[2]*y*y + params[3]*z*z) ;
    double zderiv = params[0]*2*z*params[3]*exp(params[1]*x + params[2]*y*y + params[3]*z*z) ;

    Eigen::Vector3d ans;
    //ans << yderiv, xderiv, zderiv;
    ans << xderiv, yderiv, zderiv;
    return ans;
  }
};

struct nonbonding_deriv_ModelFn{
  Eigen::Vector3d  operator() (const PerovData & data, std::vector<double> params, std::vector<double> feats){
    double x = feats[0];
    double y = feats[1];
    double z = feats[2];
    double xderiv = params[0]*y*params[1]*exp(params[1]*x + params[2]*z*z)  ;
    double yderiv = params[0]*exp(params[1]*x + params[2]*z*z) ;
    double zderiv = params[0]*y*2*z*params[2]*exp(params[1]*x + params[2]*z*z) ;

    Eigen::Vector3d ans;
    ans << xderiv, yderiv, zderiv;
    return ans;
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

    auto atom1 = data.orblut.lookup_col<Atomehr>(Orb(orb1));
    auto atom2 = data.orblut.lookup_col<Atomehr>(Orb(orb2));

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
void mpi_insert_hop_amplitudes(std::vector<double> & thop_amplitudes, HM hm, std::vector<std::vector<int>> hops, std::vector<double> params, std::vector<Eigen::Vector3d> coord){

  //std::set<Hop_idx> hopidxs=hoplut.lookup_idx(Hop(hoptype));
  //int nhops=hopidxs.size();
  int nhops = hops.size();
  if(icounthop<=ie){
    //for (auto hop_idx : hoplut.lookup_idx(Hop(hoptype))){
    for (auto hop : hops){
      if(icounthop>=is){
        int iproc=0;
        //auto orb1 = hoplut.lookup_col<Orb>(Hop_idx(hop_idx));
        //auto orb2 = hoplut.lookup_col<Orb1>(Hop_idx(hop_idx));
        //auto dx = hoplut.lookup_col<Dx>(Hop_idx(hop_idx));
        //auto dy = hoplut.lookup_col<Dy>(Hop_idx(hop_idx));
        //auto dz = hoplut.lookup_col<Dz>(Hop_idx(hop_idx));

        //std::vector<int> hop {int(dx.value), int(dy.value), int(dz.value), int(orb1.value), int(orb2.value)};

//double etime=MPI_Wtime()-stime;
//if(me==iproc) cout <<"timing hop="<<etime<<endl;

//stime=MPI_Wtime();
        complex<double> hval = hm(hop, params, coord);
//etime=MPI_Wtime()-stime;
//if(me==iproc) cout <<"timing hval="<<etime<<endl;
//stime=MPI_Wtime();
        thop_amplitudes[jcounthop]=(double) 1.0*hop[0];
        thop_amplitudes[jcounthop+1]=(double) 1.0*hop[1];
        thop_amplitudes[jcounthop+2]=(double) 1.0*hop[2];
        thop_amplitudes[jcounthop+3]=(double) 1.0*hop[3];
        thop_amplitudes[jcounthop+4]=(double) 1.0*hop[4];
        thop_amplitudes[jcounthop+5]=hval.real();
        thop_amplitudes[jcounthop+6]=hval.imag();
        jcounthop=jcounthop+7;
//double etime=MPI_Wtime()-stime;
//if(me==iproc) cout <<"timing store="<<etime<<endl;
      }
      icounthop+=1;
      if(icounthop>ie)break;
    }
  }
  else if(icounthop<is){
    icounthop+=nhops;
  }

};
/*
template <class HM>
void insert_hop_amplitudes(std::vector<std::vector<double>> & hop_amplitudes, HM hm, Hoplut hoplut, Hop hoptype, std::vector<double> params, std::vector<Eigen::Vector3d> coord){
  for (auto hop_idx : hoplut.lookup_idx(Hop(hoptype))){
    auto orb1 = hoplut.lookup_col<Orb>(Hop_idx(hop_idx));
    auto orb2 = hoplut.lookup_col<Orb1>(Hop_idx(hop_idx));
    auto dx = hoplut.lookup_col<Dx>(Hop_idx(hop_idx));
    auto dy = hoplut.lookup_col<Dy>(Hop_idx(hop_idx));
    auto dz = hoplut.lookup_col<Dz>(Hop_idx(hop_idx));

    std::vector<int> hop {int(dx.value), int(dy.value), int(dz.value), int(orb1.value), int(orb2.value)};

    complex<double> hval = hm(hop, params, coord);
    std::vector<double> hop_amplitude {1.0*dx.value, 1.0*dy.value, 1.0*dz.value, 1.0*orb1.value, 1.0*orb2.value, hval.real(), hval.imag()};
    hop_amplitudes.push_back(hop_amplitude);
  }
};
*/
void sizehopnum(std::vector<std::vector<int>> hoplut,std::vector<double> params){
//void sizehopnum(Hoplut hoplut, Hop hoptype,std::vector<double> params){
   //std::set<Hop_idx> hopidxs=hoplut.lookup_idx(Hop(hoptype)); 
   //numhop+=hopidxs.size();
   numhop+=hoplut.size();
//cout<<"rank="<<me<<" "<<hopidxs.size()<<" a "<<numhop<<endl;
};

void hopamp4mpi(){
  int nsize=numhop/nprocs-1;
  is=0;
  ndisplhoparam[0]=0;
  for (int iprocs=0;iprocs<me+1;++iprocs){
      if(iprocs!=0)is=ie+1;
      ie=is+nsize;
      if(iprocs<numhop%nprocs) ie+=1;
  }

  for (int iprocs=0;iprocs<nprocs;++iprocs){
      int nsize=numhop/nprocs;
      if(iprocs<numhop%nprocs) nsize=nsize+1;
      ndatahopamp[iprocs]=7*nsize;
      if (iprocs!=0) ndisplhoparam[iprocs]=ndatahopamp[iprocs-1]+ndisplhoparam[iprocs-1];
  }
//cout<<"pupi rank="<<me<<" "<<" is="<<is<<" ie="<<ie<<" ndata "<<ndatahopamp[me]<<" displ "<<ndisplhoparam[me]<<endl;
};
#endif
