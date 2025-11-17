#ifndef AUXPARAM_TB
#define AUXPARAM_TB


//#include "perovskite_tb.cpp"
#include "init_tb.cpp"
#include "blot/lookuptable.hpp"




  PerovData dataa=init();

  //lookuptable<Hop_idx,Hop,Orb,Orb1,Dx,Dy,Dz> hoplut;
  //int norb = (int)dataa.orblut.row_lut.size();
  int norb = dataa.norb;

  std::vector<double> pbp_params = loadtxtvec<double>("data_to_load/pbp_params.dat");
  std::vector<double> pbs_params = loadtxtvec<double>("data_to_load/pbs_params.dat");
  std::vector<double> brp_par_params = loadtxtvec<double>("data_to_load/brp_par_params.dat");
  std::vector<double> brp_perp_params = loadtxtvec<double>("data_to_load/brp_perp_params.dat");
  std::vector<double> spsig_params = loadtxtvec<double>("data_to_load/spsig_params.dat");
  std::vector<double> ppsig_params = loadtxtvec<double>("data_to_load/ppsig_params.dat");
  std::vector<double> pppi_params = loadtxtvec<double>("data_to_load/pppi_params.dat");
  std::vector<double> ppn_br_to_pb_params = loadtxtvec<double>("data_to_load/ppn_br_to_pb_params.dat");
  std::vector<double> ppn_pb_to_br_params = loadtxtvec<double>("data_to_load/ppn_pb_to_br_params.dat");
  std::vector<double> spn_pb_to_br_params = loadtxtvec<double>("data_to_load/spn_pb_to_br_params.dat");
  std::vector<double> socpb_params = loadtxtvec<double>("data_to_load/socpb_params.dat");
  std::vector<double> socbr_par_params = loadtxtvec<double>("data_to_load/socbr_par_params.dat");
  std::vector<double> socbr_perp_params = loadtxtvec<double>("data_to_load/socbr_perp_params.dat");

  HoppingModel<PerovData, pbp_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> pbp_HoppingModel(dataa);
  HoppingModel<PerovData, pbs_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> pbs_HoppingModel(dataa);
  HoppingModel<PerovData, brp_par_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> brp_par_HoppingModel(dataa);
  HoppingModel<PerovData, brp_perp_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, id_PhaseFn> brp_perp_HoppingModel(dataa);
  HoppingModel<PerovData, pppi_FrameFn, xyz_FeatureFn, pppi_ModelFn, id_PhaseFn> pppi_HoppingModel(dataa);
  HoppingModel<PerovData, ppsig_FrameFn, xyz_FeatureFn, ppsig_ModelFn, id_PhaseFn> ppsig_HoppingModel(dataa);
  HoppingModel<PerovData, spsig_FrameFn, xyz_FeatureFn, spsig_ModelFn, spsig_PhaseFn> spsig_HoppingModel(dataa);
  HoppingModel<PerovData, ppn_br_to_pb_FrameFn, xyz_FeatureFn, nonbonding_ModelFn, id_PhaseFn> ppn_br_to_pb_HoppingModel(dataa);
  HoppingModel<PerovData, ppn_pb_to_br_FrameFn, xyz_FeatureFn, nonbonding_ModelFn, id_PhaseFn> ppn_pb_to_br_HoppingModel(dataa);
  HoppingModel<PerovData, spn_pb_to_br_FrameFn, xyz_FeatureFn, nonbonding_ModelFn, id_PhaseFn> spn_pb_to_br_HoppingModel(dataa);
  HoppingModel<PerovData, socpb_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, soc_PhaseFn> socpb_HoppingModel(dataa);
  HoppingModel<PerovData, socbr_par_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, soc_PhaseFn> socbr_par_HoppingModel(dataa);
  HoppingModel<PerovData, socbr_perp_FrameFn, vcoul_FeatureFn, vcoul_ModelFn, soc_PhaseFn> socbr_perp_HoppingModel(dataa);


  double lx0;
  double ly0;
  double lz0;

  std::vector<double> latv(3);
  std::vector<double> latv0(3);

  int *ndatacoord;
  int *displcoord;
  int indexdisplcoord;


#endif
