#include "perovskite_tb.cpp"
#include "init_tb.cpp"
#include "init_hoppings.cpp"
#include "blot/lookuptable.hpp"
#include "blot/blotio.hpp"
#include "blot/tight_binding.hpp"
#include "blot/kspace.hpp"



int main(){
  PerovData data = init();
  cout << " init done" << endl;
  Hoplut hoplut(data);

  double lx = data.cell(0,0);
  double ly = data.cell(1,1);
  double lz = data.cell(2,2);
  double lx0 = lx/data.nx;
  double ly0 = ly/data.ny;
  double lz0 = lz/data.nz;

  //load parameters
  std::vector<double> pbp_params = loadtxtvec<double>("../example3/pbp_params.dat");
  std::vector<double> pbs_params = loadtxtvec<double>("../example3/pbs_params.dat");
  std::vector<double> brp_par_params = loadtxtvec<double>("../example3/brp_par_params.dat");
  std::vector<double> brp_perp_params = loadtxtvec<double>("../example3/brp_perp_params.dat");
  std::vector<double> spsig_params = loadtxtvec<double>("../example3/spsig_params.dat");
  std::vector<double> ppsig_params = loadtxtvec<double>("../example3/ppsig_params.dat");
  std::vector<double> pppi_params = loadtxtvec<double>("../example3/pppi_params.dat");
  std::vector<double> ppn_br_to_pb_params = loadtxtvec<double>("../example3/ppn_br_to_pb_params.dat");
  std::vector<double> ppn_pb_to_br_params = loadtxtvec<double>("../example3/ppn_pb_to_br_params.dat");
  std::vector<double> spn_pb_to_br_params = loadtxtvec<double>("../example3/spn_pb_to_br_params.dat");
  std::vector<double> socpb_params = loadtxtvec<double>("../example3/socpb_params.dat");
  std::vector<double> socbr_par_params = loadtxtvec<double>("../example3/socbr_par_params.dat");
  std::vector<double> socbr_perp_params = loadtxtvec<double>("../example3/socbr_perp_params.dat");

  cout << "load params done" << endl;

  //atomic coordinates
  std::vector<std::vector<double>> coord_vv = loadtxt<double>("../example3/wanndata/100fsmd/atompos.dat");
  std::vector<Eigen::Vector3d> coord;
  for (auto ii=0; ii<data.nn; ++ii){
    auto c_vv = coord_vv[ii];
    Eigen::Vector3d c0 { c_vv[0]*lx, c_vv[1]*ly, c_vv[2]*lz};
    Eigen::Vector3d c1 = translate_posvec0(data.atomlut, Atom(ii), c0, lx0, ly0, lz0);

    coord.push_back(c1);
  }


  cout << "atom coords done" << endl;

  //calculate hopping values
  
  HoppingModel<PerovData, pbp_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> pbp_HoppingModel(data);
  HoppingModel<PerovData, pbs_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> pbs_HoppingModel(data);
  HoppingModel<PerovData, brp_par_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> brp_par_HoppingModel(data);
  HoppingModel<PerovData, brp_perp_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> brp_perp_HoppingModel(data);
  HoppingModel<PerovData, pppi_FrameFn, xyz_FeatureFn, pppi_ModelFn, id_PhaseFn> pppi_HoppingModel(data);
  HoppingModel<PerovData, ppsig_FrameFn, xyz_FeatureFn, ppsig_ModelFn, id_PhaseFn> ppsig_HoppingModel(data);
  HoppingModel<PerovData, spsig_FrameFn, xyz_FeatureFn, spsig_ModelFn, spsig_PhaseFn> spsig_HoppingModel(data);
  HoppingModel<PerovData, ppn_br_to_pb_FrameFn, xyz_FeatureFn, nonbonding_ModelFn, id_PhaseFn> ppn_br_to_pb_HoppingModel(data);
  HoppingModel<PerovData, ppn_pb_to_br_FrameFn, xyz_FeatureFn, nonbonding_ModelFn, id_PhaseFn> ppn_pb_to_br_HoppingModel(data);
  HoppingModel<PerovData, spn_pb_to_br_FrameFn, xyz_FeatureFn, nonbonding_ModelFn, id_PhaseFn> spn_pb_to_br_HoppingModel(data);
  HoppingModel<PerovData, socpb_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, soc_PhaseFn> socpb_HoppingModel(data);
  HoppingModel<PerovData, socbr_par_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, soc_PhaseFn> socbr_par_HoppingModel(data);
  HoppingModel<PerovData, socbr_perp_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, soc_PhaseFn> socbr_perp_HoppingModel(data);

  cout << "set up hopping models done" << endl;
  std::vector<std::vector<double>> hop_amplitudes;

  insert_hop_amplitudes(hop_amplitudes, pbp_HoppingModel, hoplut.pbp, pbp_params, coord);
  insert_hop_amplitudes(hop_amplitudes, pbs_HoppingModel, hoplut.pbs, pbs_params, coord);
  insert_hop_amplitudes(hop_amplitudes, brp_par_HoppingModel, hoplut.brp_par, brp_par_params, coord);
  insert_hop_amplitudes(hop_amplitudes, brp_perp_HoppingModel, hoplut.brp_perp, brp_perp_params, coord);
  insert_hop_amplitudes(hop_amplitudes, pppi_HoppingModel, hoplut.pppi, pppi_params, coord);
  insert_hop_amplitudes(hop_amplitudes, ppsig_HoppingModel, hoplut.ppsig, ppsig_params, coord);
  insert_hop_amplitudes(hop_amplitudes, spsig_HoppingModel, hoplut.spsig, spsig_params, coord);
  insert_hop_amplitudes(hop_amplitudes, ppn_br_to_pb_HoppingModel, hoplut.ppn_br_to_pb, ppn_br_to_pb_params, coord);
  insert_hop_amplitudes(hop_amplitudes, ppn_pb_to_br_HoppingModel, hoplut.ppn_pb_to_br, ppn_pb_to_br_params, coord);
  insert_hop_amplitudes(hop_amplitudes, spn_pb_to_br_HoppingModel, hoplut.spn_pb_to_br, spn_pb_to_br_params, coord);
  insert_hop_amplitudes(hop_amplitudes, socpb_HoppingModel, hoplut.socpb, socpb_params, coord);
  insert_hop_amplitudes(hop_amplitudes, socbr_par_HoppingModel, hoplut.socbr_par, socbr_par_params, coord);
  insert_hop_amplitudes(hop_amplitudes, socbr_perp_HoppingModel, hoplut.socbr_perp, socbr_perp_params, coord);

  savetxt<double>("hopping.dat",hop_amplitudes);

  cout << "hopping amplitudes done" << endl;

  //calculate energy bands
  SlepcInitWrap();
  //set k-path
  vector<vector<double> > kp0 = { {0,0,0}, {0.5,0,0}, {0.5,0.5,0}};
  auto kp = kpath(kp0, 10);

  TBHammaker hm = {3, 208, hop_amplitudes};
  cout << "kpath done" << endl;

  //set options
  //Solve_ham_interval_options op = {-10.0,10.0};
  //Solve_ham_interval_options op = {-3.0,3.0};
  //Solve_ham_target_options op = {0.0,8};
  Solve_ham_dense_options op;
  //

  cout << "options done" << endl;
  //solve
  solve_ham_params(hm, kp, false, op);


  return 0;
}
