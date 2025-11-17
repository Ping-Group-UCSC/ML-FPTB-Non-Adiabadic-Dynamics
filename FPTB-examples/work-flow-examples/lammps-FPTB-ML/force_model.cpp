#include "perovskite_tb.cpp"
#include <math.h>
#include <cmath>
#include "blot/blotmath.hpp"
#include "blot/units.hpp"

//template for a functor that computes force vectors given atomic positions and features
//

template<class SystemData, class ModelFn, class PhaseFn >
struct ForceModel{
  SystemData & data;
  ModelFn derivmodelfn;
  PhaseFn phasefn;

  ForceModel(SystemData & data1) : data(data1) { };

  Eigen::Vector3cd operator() (std::vector<int> hop,  std::vector<double> params, std::vector<double> feats, const Eigen::MatrixXd & dfeats,  std::vector<Eigen::Vector3d> atompos){
    auto fval = derivmodelfn(data,params,feats);
    Eigen::Vector3cd fval1 = (dfeats * fval) /phasefn(data,hop);

    return fval1;
  }
};

void set_force_cart(Eigen::Vector3cd force_value, std::vector<int> hop, MatWrap& fx, MatWrap& fy, MatWrap& fz, int row_start_x, int row_end_x, int row_start_y, int row_end_y, int row_start_z, int row_end_z){

  int from = hop[3];
  int to = hop[4];

    MatSetValue(fx.A, from-1, to-1, force_value(0) , ADD_VALUES);
    if (from!=to) MatSetValue(fx.A, to-1, from-1, conj(force_value(0)) , ADD_VALUES);

    MatSetValue(fy.A, from-1, to-1, force_value(1) , ADD_VALUES);
    if (from!=to) MatSetValue(fy.A, to-1, from-1, conj(force_value(1)) , ADD_VALUES);

    MatSetValue(fz.A, from-1, to-1, force_value(2) , ADD_VALUES);
    if (from!=to) MatSetValue(fz.A, to-1, from-1, conj(force_value(2)) , ADD_VALUES);
}



//calculate hellmann-feynman forces 
// [atompos] = angstrom
// [forces] = N 
std::vector<Eigen::Vector3d> force_model(PerovData& data, std::vector<Eigen::Vector3d>& atompos, VecWrap& psi0, std::vector<int>& iorbs){

  vcoul_cutoff_DFeatureFn vcoul_cutoff_dfeatfn;
  xyz_DFeatureFn xyz_dfeatfn;
  xyz_FeatureFn xyz_featfn;

  pbs_FrameFn pbs_framefn;
  spsig_FrameFn spsig_framefn;
  ppsig_FrameFn ppsig_framefn;
  pppi_FrameFn pppi_framefn;
  ppn_br_to_pb_FrameFn ppn_br_to_pb_framefn;
  ppn_pb_to_br_FrameFn ppn_pb_to_br_framefn;
  spn_pb_to_br_FrameFn spn_pb_to_br_framefn;

  ForceModel<PerovData, vcoul_deriv_ModelFn, id_PhaseFn> fcoul_ForceModel(data);
  ForceModel<PerovData, vcoul_deriv_ModelFn, soc_PhaseFn> soc_ForceModel(data);
  ForceModel<PerovData, spsig_deriv_ModelFn, spsig_PhaseFn> spsig_ForceModel(data);
  ForceModel<PerovData, ppsig_deriv_ModelFn, id_PhaseFn> ppsig_ForceModel(data);
  ForceModel<PerovData, pppi_deriv_ModelFn, id_PhaseFn> pppi_ForceModel(data);
  ForceModel<PerovData, nonbonding_deriv_ModelFn, id_PhaseFn> ppn_br_to_pb_ForceModel(data);
  ForceModel<PerovData, nonbonding_deriv_ModelFn, id_PhaseFn> ppn_pb_to_br_ForceModel(data);
  ForceModel<PerovData, nonbonding_deriv_ModelFn, id_PhaseFn> spn_pb_to_br_ForceModel(data);

  std::vector<double> tans_x(ndatacoord[me]);
  std::vector<double> ans_temp(3*data.nn);

  VecWrap psi = psi0.get(iorbs);

  int iatom;
  Eigen::Vector3d temp_force(3);

  for (int i=0; i<ndatacoord[me]/3;++i){

    iatom = i+(displcoord[me]/3);
    Atomehr atom1(iatom);

    auto ix = data.atomlut.lookup_col<Atomehr_ix>(atom1);
    auto iy = data.atomlut.lookup_col<Atomehr_iy>(atom1);
    auto iz = data.atomlut.lookup_col<Atomehr_iz>(atom1);
    auto atom1_typ = data.atomlut.lookup_col<Atomehr_typ>(atom1);
  
    MatWrap fx(data.norb,data.norb,true);
    MatWrap fy(data.norb,data.norb,true);
    MatWrap fz(data.norb,data.norb,true);
  
    int row_start_x, row_end_x;
    MatGetOwnershipRange(fx.A, &row_start_x, &row_end_x);
    int row_start_y, row_end_y;
    MatGetOwnershipRange(fy.A, &row_start_y, &row_end_y);
    int row_start_z, row_end_z;
    MatGetOwnershipRange(fz.A, &row_start_z, &row_end_z);

    //derivative of onsite and soc terms
    //add deriv of atom1
    std::vector<std::vector<int>> nbrs = data.vcutlut[atom1.value];
    std::vector<int> selfnbr {atom1.value, atom1.value, 0, 0, 0};
    nbrs.push_back(selfnbr);
    for (auto nbr : nbrs){
      auto dx = nbr[2];
      auto dy = nbr[3];
      auto dz = nbr[4];
      auto ia2 = nbr[1];

      auto atom2 = Atomehr(ia2);
      auto atom_typ2 = data.atomlut.lookup_col<Atomehr_typ>(Atomehr(atom2));

      //calculate fcoul features
      Eigen::MatrixXd vcoul_cutoff_dfeats = vcoul_cutoff_dfeatfn(data, atompos, atom2, atom1, Dx(dx), Dy(dy), Dz(dz));
      
      //derivative of onsite terms
      for (auto orb2 : data.orblut.lookup_idx(Atomehr(atom2))){
        auto orb_typ = data.orblut.lookup_col<Orb_typ>(Orb(orb2));
        double scalefac = 0.0;
        if (orb_typ == "Pb:p") scalefac = data.pbp_params[0];
        if (orb_typ == "Pb:s") scalefac = data.pbs_params[0];
        if (orb_typ == "Br:p"){
          auto ddir = data.orblut.lookup_col<Ddir>(orb2);
          auto pdir = data.orblut.lookup_col<Pdir>(orb2);
          if (ddir==pdir){
            scalefac = data.brp_par_params[0];
          }
          else{
            scalefac = data.brp_perp_params[0];
          }
        }
    
        std::vector<double> params = {scalefac};
        std::vector<double> feats = {} ;
        std::vector<int> hop = { 0, 0, 0, orb2.value, orb2.value};
        Eigen::Vector3cd vcoul_deriv_value = fcoul_ForceModel(hop, params, feats, vcoul_cutoff_dfeats , atompos);

        set_force_cart(vcoul_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

      }


      //derivative of soc terms
      if(atom_typ2=="Pb"){
        int orb1 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("x") );
        int orb2 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("x") );
        int orb3 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("y") );
        int orb4 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("y") );
        int orb5 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("z") );
        int orb6 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("z") );

        std::vector<int> h1 = { 0, 0, 0, orb1, orb3}; 
        std::vector<int> h2 = { 0, 0, 0, orb2, orb4}; 
        std::vector<int> h3 = { 0, 0, 0, orb2, orb5}; 
        std::vector<int> h4 = { 0, 0, 0, orb1, orb6}; 
        std::vector<int> h5 = { 0, 0, 0, orb3, orb6}; 
        std::vector<int> h6 = { 0, 0, 0, orb4, orb5}; 

        double scalefac = data.socpb_params[0];

        std::vector<double> params = {scalefac};
        std::vector<double> feats = {} ;

        Eigen::Vector3cd socpb_deriv_value_h1= soc_ForceModel(h1, params, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h2= soc_ForceModel(h2, params, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h3= soc_ForceModel(h3, params, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h4= soc_ForceModel(h4, params, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h5= soc_ForceModel(h5, params, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h6= soc_ForceModel(h6, params, feats, vcoul_cutoff_dfeats, atompos);

        set_force_cart(socpb_deriv_value_h1, h1, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h2, h2, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h3, h3, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h4, h4, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h5, h5, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h6, h6, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
      }
        
      if(atom_typ2=="Br1"){
        double scalefac_par = data.socbr_par_params[0];
        double scalefac_perp = data.socbr_perp_params[0];
        std::vector<double> params_par = {scalefac_par};
        std::vector<double> params_perp = {scalefac_perp};

        int orb1 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("x") );
        int orb2 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("x") );
        int orb3 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("y") );
        int orb4 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("y") );
        int orb5 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("z") );
        int orb6 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("z") );

        std::vector<int> h1 = { 0, 0, 0, orb1, orb3}; 
        std::vector<int> h2 = { 0, 0, 0, orb2, orb4}; 
        std::vector<int> h3 = { 0, 0, 0, orb2, orb5}; 
        std::vector<int> h4 = { 0, 0, 0, orb1, orb6}; 
        std::vector<int> h5 = { 0, 0, 0, orb3, orb6}; 
        std::vector<int> h6 = { 0, 0, 0, orb4, orb5}; 

        std::vector<double> feats = {} ;

        Eigen::Vector3cd socpb_deriv_value_h1= soc_ForceModel(h1, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h2= soc_ForceModel(h2, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h3= soc_ForceModel(h3, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h4= soc_ForceModel(h4, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h5= soc_ForceModel(h5, params_perp, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h6= soc_ForceModel(h6, params_perp, feats, vcoul_cutoff_dfeats, atompos);

        set_force_cart(socpb_deriv_value_h1, h1, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h2, h2, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h3, h3, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h4, h4, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h5, h5, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h6, h6, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
      }

      if(atom_typ2=="Br2"){
        double scalefac_par = data.socbr_par_params[0];
        double scalefac_perp = data.socbr_perp_params[0];
        std::vector<double> params_par = {scalefac_par};
        std::vector<double> params_perp = {scalefac_perp};

        int orb1 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("x") );
        int orb2 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("x") );
        int orb3 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("y") );
        int orb4 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("y") );
        int orb5 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("z") );
        int orb6 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("z") );

        std::vector<int> h1 = { 0, 0, 0, orb1, orb3}; 
        std::vector<int> h2 = { 0, 0, 0, orb2, orb4}; 
        std::vector<int> h3 = { 0, 0, 0, orb2, orb5}; 
        std::vector<int> h4 = { 0, 0, 0, orb1, orb6}; 
        std::vector<int> h5 = { 0, 0, 0, orb3, orb6}; 
        std::vector<int> h6 = { 0, 0, 0, orb4, orb5}; 

        std::vector<double> feats = {} ;

        Eigen::Vector3cd socpb_deriv_value_h1= soc_ForceModel(h1, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h2= soc_ForceModel(h2, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h3= soc_ForceModel(h3, params_perp, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h4= soc_ForceModel(h4, params_perp, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h5= soc_ForceModel(h5, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h6= soc_ForceModel(h6, params_par, feats, vcoul_cutoff_dfeats, atompos);

        set_force_cart(socpb_deriv_value_h1, h1, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h2, h2, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h3, h3, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h4, h4, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h5, h5, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h6, h6, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
      }

      if(atom_typ2=="Br3"){
        double scalefac_par = data.socbr_par_params[0];
        double scalefac_perp = data.socbr_perp_params[0];
        std::vector<double> params_par = {scalefac_par};
        std::vector<double> params_perp = {scalefac_perp};

        int orb1 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("x") );
        int orb2 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("x") );
        int orb3 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("y") );
        int orb4 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("y") );
        int orb5 = data.pbp_idx(Atomehr(atom2), Spin("up"), Pdir("z") );
        int orb6 = data.pbp_idx(Atomehr(atom2), Spin("dn"), Pdir("z") );

        std::vector<int> h1 = { 0, 0, 0, orb1, orb3}; 
        std::vector<int> h2 = { 0, 0, 0, orb2, orb4}; 
        std::vector<int> h3 = { 0, 0, 0, orb2, orb5}; 
        std::vector<int> h4 = { 0, 0, 0, orb1, orb6}; 
        std::vector<int> h5 = { 0, 0, 0, orb3, orb6}; 
        std::vector<int> h6 = { 0, 0, 0, orb4, orb5}; 

        std::vector<double> feats = {} ;

        Eigen::Vector3cd socpb_deriv_value_h1= soc_ForceModel(h1, params_perp, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h2= soc_ForceModel(h2, params_perp, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h3= soc_ForceModel(h3, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h4= soc_ForceModel(h4, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h5= soc_ForceModel(h5, params_par, feats, vcoul_cutoff_dfeats, atompos);
        Eigen::Vector3cd socpb_deriv_value_h6= soc_ForceModel(h6, params_par, feats, vcoul_cutoff_dfeats, atompos);

        set_force_cart(socpb_deriv_value_h1, h1, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h2, h2, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h3, h3, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h4, h4, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h5, h5, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
        set_force_cart(socpb_deriv_value_h6, h6, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);
      }

    }

    //derivative of nearest neighbor hopping terms
    //

    Atomehr atomPb;
    Atomehr atomBr;
    Atomehr atom2;
    std::set<Shell_idx> nnshells;
    if (atom1_typ == "Pb"){ 
      nnshells = data.shelllut.lookup_idx(Atomehr(atom1.value),Shell("pbbr_shell_1")) ;
    }
    else if(atom1_typ =="Br"){
      nnshells = data.shelllut.lookup_idx(Atomehr_nbr(atom1.value),Shell("pbbr_shell_1")) ;
    }
    
    for(auto ii : nnshells){
      auto atompb = data.shelllut.lookup_col<Atomehr>(Shell_idx(ii));
      auto atombr_nbr = data.shelllut.lookup_col<Atomehr_nbr>(Shell_idx(ii));
      auto atombr = Atomehr(atombr_nbr.value);

      if (atom1_typ == "Pb"){ 
        atom2 = atombr;
      }
      else if (atom1_typ == "Br"){
        atom2 = atompb;
      }

      auto dx = data.shelllut.lookup_col<Dx>(Shell_idx(ii));
      auto dy = data.shelllut.lookup_col<Dy>(Shell_idx(ii));
      auto dz = data.shelllut.lookup_col<Dz>(Shell_idx(ii));

      for( auto orbpb : data.orblut.lookup_idx(Atomehr(atompb))){
        for( auto orbbr : data.orblut.lookup_idx(Atomehr(atombr.value))){

          auto ddirbr = data.orblut.lookup_col<Ddir>(orbbr);

          auto pdirpb = data.orblut.lookup_col<Pdir>(orbpb);
          auto pdirbr = data.orblut.lookup_col<Pdir>(orbbr);

          auto spinpb = data.orblut.lookup_col<Spin>(orbpb);
          auto spinbr = data.orblut.lookup_col<Spin>(orbbr);
          if (data.orblut.lookup_col<Orb_typ>(orbpb) == "Pb:s"){
            if (ddirbr==pdirbr && spinpb==spinbr){
              std::vector<int> hop = {dx.value, dy.value, dz.value, orbpb.value, orbbr.value};

              Eigen::MatrixXd frame = spsig_framefn(data, hop, atompos);
              std::vector<double> xyz_feat = xyz_featfn(data, frame, hop, atompos);

              Eigen::MatrixXd dfeat = xyz_dfeatfn(data, frame, hop, atompos, atom2, Dx(0), Dy(0), Dz(0));

              Eigen::Vector3cd spsig_deriv_value = spsig_ForceModel(hop, data.spsig_params, xyz_feat, dfeat, atompos);
              set_force_cart(spsig_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

            }
            else if(spinpb==spinbr){
              std::vector<int> hop = {dx.value, dy.value, dz.value, orbpb.value, orbbr.value};
              Eigen::MatrixXd frame = spn_pb_to_br_framefn(data, hop, atompos);
              std::vector<double> xyz_feat = xyz_featfn(data, frame, hop, atompos);

              Eigen::MatrixXd dfeat = xyz_dfeatfn(data, frame, hop, atompos, atom2, Dx(0), Dy(0), Dz(0));
              
              Eigen::Vector3cd spn_pb_to_br_deriv_value = spn_pb_to_br_ForceModel(hop, data.spn_pb_to_br_params, xyz_feat, dfeat, atompos);
              set_force_cart(spn_pb_to_br_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

            }
          }
          else if (data.orblut.lookup_col<Orb_typ>(orbpb) == "Pb:p"){
            if (ddirbr==pdirpb && spinpb==spinbr){
              if (ddirbr==pdirbr){
                std::vector<int> hop = {dx.value, dy.value, dz.value, orbpb.value, orbbr.value};
                Eigen::MatrixXd frame = ppsig_framefn(data, hop, atompos);
                std::vector<double> xyz_feat = xyz_featfn(data, frame, hop, atompos);

                Eigen::MatrixXd dfeat = xyz_dfeatfn(data, frame, hop, atompos, atom2, Dx(0), Dy(0), Dz(0));
                
                Eigen::Vector3cd ppsig_deriv_value = ppsig_ForceModel(hop, data.ppsig_params, xyz_feat, dfeat, atompos);
                set_force_cart(ppsig_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

              }
              else{
                std::vector<int> hop = {dx.value, dy.value, dz.value, orbpb.value, orbbr.value};
                Eigen::MatrixXd frame = ppn_pb_to_br_framefn(data, hop, atompos);
                std::vector<double> xyz_feat = xyz_featfn(data, frame, hop, atompos);

                Eigen::MatrixXd dfeat = xyz_dfeatfn(data, frame, hop, atompos, atom2, Dx(0), Dy(0), Dz(0));

                Eigen::Vector3cd ppn_pb_to_br_deriv_value = ppn_pb_to_br_ForceModel(hop, data.ppn_pb_to_br_params, xyz_feat, dfeat, atompos);
                set_force_cart(ppn_pb_to_br_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

              }
            }
            else if (ddirbr!=pdirpb && spinpb==spinbr){
              if (pdirpb==pdirbr){
                std::vector<int> hop = {dx.value, dy.value, dz.value, orbpb.value, orbbr.value};

                Eigen::MatrixXd frame = pppi_framefn(data, hop, atompos);
                std::vector<double> xyz_feat = xyz_featfn(data, frame, hop, atompos);

                Eigen::MatrixXd dfeat = xyz_dfeatfn(data, frame, hop, atompos, atom2, Dx(0), Dy(0), Dz(0));

                Eigen::Vector3cd pppi_deriv_value = pppi_ForceModel(hop, data.pppi_params, xyz_feat, dfeat, atompos);
                set_force_cart(pppi_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);

              }
              else if (ddirbr==pdirbr){
                std::vector<int> hop = {-dx.value, -dy.value, -dz.value, orbbr.value, orbpb.value};
                Eigen::MatrixXd frame = ppn_br_to_pb_framefn(data, hop, atompos);
                std::vector<double> xyz_feat = xyz_featfn(data, frame, hop, atompos);

                Eigen::MatrixXd dfeat = xyz_dfeatfn(data, frame, hop, atompos, atom2, Dx(0), Dy(0), Dz(0));

                Eigen::Vector3cd ppn_br_to_pb_deriv_value = ppn_br_to_pb_ForceModel(hop, data.ppn_br_to_pb_params, xyz_feat, dfeat, atompos);
                set_force_cart(ppn_br_to_pb_deriv_value, hop, fx, fy, fz, row_start_x, row_end_x, row_start_y, row_end_y, row_start_z, row_end_z);


              }
            }
          }
        }
      }
    }

    //force expectation value
    MatAssemblyBegin(fx.A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fy.A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(fz.A, MAT_FINAL_ASSEMBLY);

    MatAssemblyEnd(fx.A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fy.A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(fz.A, MAT_FINAL_ASSEMBLY);

    VecWrap fxpsi(PETSC_COMM_SELF, psi.n);
    VecWrap fypsi(PETSC_COMM_SELF, psi.n);
    VecWrap fzpsi(PETSC_COMM_SELF, psi.n);

    PetscScalar fx_val, fy_val, fz_val;

    MatMult(fx.A, psi.V, fxpsi.V);
    MatMult(fy.A, psi.V, fypsi.V);
    MatMult(fz.A, psi.V, fzpsi.V);

    VecDot(psi.V, fxpsi.V, &fx_val);
    VecDot(psi.V, fypsi.V, &fy_val);
    VecDot(psi.V, fzpsi.V, &fz_val);

    temp_force << real(fx_val)*elemSI/angSI , real(fy_val)*elemSI/angSI , real(fz_val)*elemSI/angSI ;

    tans_x[3*i] = temp_force(0);
    tans_x[(3*i)+1] = temp_force(1);
    tans_x[(3*i)+2] = temp_force(2);
  }


  MPI_Allgatherv(tans_x.data(),ndatacoord[me],MPI_DOUBLE,ans_temp.data(),ndatacoord,displcoord, MPI_DOUBLE, MPI_COMM_WORLD);

  std::vector<Eigen::Vector3d> ans(data.nn, Eigen::Vector3d(3));
  for (int i=0; i<data.nn;++i){
    ans[i] << ans_temp[3*i], ans_temp[(3*i)+1], ans_temp[(3*i)+2];
  }
  return ans;

}

