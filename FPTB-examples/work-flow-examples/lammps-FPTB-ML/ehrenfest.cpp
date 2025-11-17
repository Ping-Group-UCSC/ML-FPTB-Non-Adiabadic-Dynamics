#include <mpi.h>
#include <math.h>
#include <fstream>
#include <regex>

#include "lammps.h"

#include "atom.h"
#include "library.h"

#include "modify.h"
#include "fix.h"
#include "fix_external.h"


#include "auxparam_tb.hpp"
#include "auxparam.hpp"

#include "perovskite_tb.cpp"
#include "force_model.cpp"
#include "diffusion.cpp"
#include "charges.cpp"
#include "init_hoppings.cpp"
#include "blot/lookuptable.hpp"
#include "blot/blotio.hpp"
#include "blot/tight_binding.hpp"
#include "blot/dynamics.hpp"
#include "blot/units.hpp"


using namespace LAMMPS_NS;

void ehrenfest_callback (void *, bigint, int, int *, double **, double **);


static char help[] = "ehrenfest example.\n\n";

Hoplut hoplut;
int main(int narg, char **arg){


  MPI_Init(&narg,&arg);

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  PetscInitialize(&narg,&arg,(char*)0,help);

  lmp = NULL;
  lmp = new LAMMPS(0,NULL,MPI_COMM_WORLD);

  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[1],"r");
    if (fp == NULL) {
      printf("ERROR: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  int n;
  char line[1024];
  while (1) {
    if (me == 0) {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
      if (n == 0) fclose(fp);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
    lammps_command(lmp,line);
  }


  hoplut.setup(dataa);

  lx0 = dataa.cell(0,0)/dataa.nx;
  ly0 = dataa.cell(1,1)/dataa.ny;
  lz0 = dataa.cell(2,2)/dataa.nz;

  latv[0]=dataa.cell(0,0);
  latv[1]=dataa.cell(1,1);
  latv[2]=dataa.cell(2,2);
  latv0[0]=lx0;
  latv0[1]=ly0;
  latv0[2]=lz0;

  ndatacoord = NULL;
  displcoord = NULL;
  ndatacoord= new int[nprocs];
  displcoord= new int[nprocs];
  ndatacd = NULL;
  displcd = NULL;
  ndatacd= new int[nprocs];
  displcd= new int[nprocs];
  displcoord[0]= 0;
  displcd[0]= 0;
  indexdisplcoord=0;
  iscd=0;
  for (int iprocs=0;iprocs<nprocs;++iprocs){
      int nsize=dataa.nn/nprocs;
      if(iprocs<dataa.nn%nprocs) nsize=nsize+1;
      if(iprocs<me) indexdisplcoord=indexdisplcoord+nsize;
      ndatacoord[iprocs]=3*nsize;
      ndatacd[iprocs]=nsize;
      if (iprocs!=0) displcoord[iprocs]=ndatacoord[iprocs-1]+displcoord[iprocs-1];
      if (iprocs!=0) displcd[iprocs]=ndatacd[iprocs-1]+displcd[iprocs-1];
      if(iprocs<me) iscd=iecd+1;
      if(iprocs<me+1) iecd=iscd+nsize;
  }


  numhop=0;
  sizehopnum(hoplut.pbp,pbp_params);
  sizehopnum(hoplut.pbs,pbs_params);
  sizehopnum(hoplut.brp_par,brp_par_params);
  sizehopnum(hoplut.brp_perp,brp_perp_params);
  sizehopnum(hoplut.pppi,pppi_params);
  sizehopnum(hoplut.ppsig,ppsig_params);
  sizehopnum(hoplut.spsig,spsig_params);
  sizehopnum(hoplut.ppn_br_to_pb,ppn_br_to_pb_params);
  sizehopnum(hoplut.ppn_pb_to_br,ppn_pb_to_br_params);
  sizehopnum(hoplut.spn_pb_to_br,spn_pb_to_br_params);
  sizehopnum(hoplut.socpb,socpb_params);
  sizehopnum(hoplut.socbr_par,socbr_par_params);
  sizehopnum(hoplut.socbr_perp,socbr_perp_params);

  ndatahopamp= NULL;
  ndisplhoparam= NULL;
  ndatahopamp= new int[nprocs];
  ndisplhoparam= new int[nprocs];
                
  hopamp4mpi();

  tewald_phi=NULL;
  tewald_phi= new double[dataa.nn];

  coefewald_field=elemSI*elemSI*(1.0/4/M_PI/eps0SI)/angSI/angSI;

  int ifix = lmp->modify->find_fix("4");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_callback(ehrenfest_callback,&lmp);


  lmp->input->one("run 20");
  lmp->input->one("write_restart nvt20.pvk.restart");



  delete[] ndatacoord;
  delete[] displcoord;
  delete[] ndatacd;
  delete[] displcd;
  delete[] ndatahopamp;
  delete[] ndisplhoparam;
  delete[] tewald_phi;


  MPI_Barrier(MPI_COMM_WORLD);
  PetscFinalize(); 
  MPI_Finalize();

  return 0;
}


void ehrenfest_callback(void *ptr, bigint ntimestep,
                    int nlocal, int *id, double **x, double **f)
{

  double *xx = NULL; 

  xx = new double[3*dataa.nn];
  lammps_gather_atoms(lmp,(char *) "x",1,3,xx);

    

  Eigen::VectorXd coordarray(dataa.nn*3);
  std::vector<double> c_vv(3);
  Eigen::VectorXd temp(ndatacoord[me]);
  for (int ii=0; ii<ndatacoord[me]/3; ++ii){
    for (int i=0;i<3;++i){
      c_vv[i]=xx[i+(3*ii+displcoord[me])];
        if (latv[i]-c_vv[i] < latv0[i]/2){
        c_vv[i] -=latv[i];
        }
    }
    Eigen::Vector3d c0 { c_vv[0], c_vv[1], c_vv[2]};
    Eigen::Vector3d c1 = translate_posvec0(dataa.atomlut, Atomehr(ii+indexdisplcoord), c0, lx0, ly0, lz0);
    for (int i=0;i<3;++i){
      temp(ii*3+i) = c1[i];
    }
  }


  MPI_Allgatherv(temp.data(),ndatacoord[me],MPI_DOUBLE,coordarray.data(),ndatacoord,displcoord, MPI_DOUBLE, MPI_COMM_WORLD);

  std::vector<Eigen::Vector3d> coord(dataa.nn);
  for (int ii=0;ii<dataa.nn; ++ii){
    coord[ii]=Eigen::Vector3d {coordarray[ii*3], coordarray[1+ii*3], coordarray[2+ii*3]};
  }
  // write positions used to calculate hopping amplitudes to file
  // NOTE file "pos-ehrn.xyz" will be appended no rewritten if run multiple times 
  if (me==0){
    std::ofstream out;
    out.open("pos-ehrn.xyz",std::ios::app);
    out << ntimestep << endl;

    for (int ii=0;ii<dataa.nn; ++ii){
      out << coordarray[ii*3] << " " << coordarray[1+ii*3] << " " << coordarray[2+ii*3] << endl;
      }
    out.close();
    }

  cell=dataa.cell;
  gvecs=2.0*M_PI*cell.inverse();
  vol = abs(cell.determinant());
  sqrt2M_PIsig=sqrt(2.0/M_PI)/sigewald; // sigewald is declared as a constant in auxparam.hpp
  fM_PIvol=4.0*M_PI/vol;
  hsigsig=-0.5*sigewald*sigewald;
  onesqrt2sig=1.0/sqrt(2.0)/sigewald;

  std::vector<double> tewald_temp(ndatacd[me]);
  int ewald_i;

  for (int iatom=displcd[me];iatom<ndatacd[me]+displcd[me];++iatom){
    ewald_i= iatom - displcd[me];
    tewald_temp[ewald_i]=coul_cutoff_AG(dataa, iatom, coord,2,2,2);
  }
  MPI_Allgatherv(tewald_temp.data(),ndatacd[me],MPI_DOUBLE,tewald_phi,ndatacd,displcd, MPI_DOUBLE, MPI_COMM_WORLD);

  icounthop=0;
  jcounthop=0;


  std::vector<std::vector<double>> hop_amplitudes;
  std::vector<double> thop_amplitudes(ndatahopamp[me]);
  std::vector<double> hoparray(7*numhop);


  mpi_insert_hop_amplitudes(thop_amplitudes, pbp_HoppingModel, hoplut.pbp, pbp_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, pbs_HoppingModel, hoplut.pbs, pbs_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, brp_par_HoppingModel, hoplut.brp_par, brp_par_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, brp_perp_HoppingModel, hoplut.brp_perp, brp_perp_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, pppi_HoppingModel, hoplut.pppi, pppi_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, ppsig_HoppingModel, hoplut.ppsig, ppsig_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, spsig_HoppingModel, hoplut.spsig, spsig_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, ppn_br_to_pb_HoppingModel, hoplut.ppn_br_to_pb, ppn_br_to_pb_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, ppn_pb_to_br_HoppingModel, hoplut.ppn_pb_to_br, ppn_pb_to_br_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, spn_pb_to_br_HoppingModel, hoplut.spn_pb_to_br, spn_pb_to_br_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, socpb_HoppingModel, hoplut.socpb, socpb_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, socbr_par_HoppingModel, hoplut.socbr_par, socbr_par_params, coord);
  mpi_insert_hop_amplitudes(thop_amplitudes, socbr_perp_HoppingModel, hoplut.socbr_perp, socbr_perp_params, coord);

  MPI_Allgatherv(thop_amplitudes.data(),ndatahopamp[me],MPI_DOUBLE,hoparray.data(),ndatahopamp,ndisplhoparam, MPI_DOUBLE, MPI_COMM_WORLD);


  for (int ii=0;ii<numhop; ++ii){
    int jj=ii*7;
    std::vector<double> c3 {hoparray[jj], hoparray[1+jj], hoparray[2+jj],hoparray[3+jj],hoparray[4+jj],hoparray[5+jj],hoparray[6+jj]};
    hop_amplitudes.push_back(c3);
  }

  if(me==0)  savetxt<double>("hopping.dat",hop_amplitudes);

  TBHammaker hm = {3, norb, hop_amplitudes};
  MatWrap h(norb,norb);
  hm.create(h, std::vector<double>(3,0.0));
  ETRS_static_stepper stepper(h, 4);

  auto const atompb1 = *dataa.atomlut.lookup_idx(Atomehr_typ("Pb"),Atomehr_ix(1), Atomehr_iy(1), Atomehr_iz(1)).begin();
  auto const atombr1 = *dataa.atomlut.lookup_idx(Atomehr_typ("Br1"),Atomehr_ix(2), Atomehr_iy(1), Atomehr_iz(1)).begin();
  auto const atompb2 = *dataa.atomlut.lookup_idx(Atomehr_typ("Pb"),Atomehr_ix(2), Atomehr_iy(1), Atomehr_iz(1)).begin();

  auto const orbpb1p = *dataa.orblut.lookup_idx(Orb_typ("Pb:p"), Atomehr(atompb1), Spin("up")).begin();
  auto const orbbr1px = *dataa.orblut.lookup_idx(Orb_typ("Br:p"), Atomehr(atombr1), Pdir("x"), Spin("up")).begin();
  auto const orbpb2s = *dataa.orblut.lookup_idx(Orb_typ("Pb:s"), Atomehr(atompb2), Spin("up")).begin();
  std::vector<int> idxseq_v = {orbpb1p.value-1, orbbr1px.value-1, orbpb2s.value-1};
  Simple_monitor monitor(1,"trajectory.dat", idxseq_v);

  
  //time step
  double dt = 1.0e-6; 
  int nt = 500;

  VecWrap psi0(norb);

  std::vector<int> iorbs ;
  for (int iorb=0; iorb<norb; ++iorb){
    iorbs.push_back(iorb);
  }
  
  // set up initial values for psi if first timestep
  if (ntimestep==0){ 
    std::vector<std::complex<double>> vecin;
    double real;
    double imag;
    if (me==0) cout << "ORIGINAL" << endl;
    std::regex rgx("^(?=[iI.\\d+-])([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?(?![iI.\\d]))?([+-]?(?:(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)?[iI])?$");
    std::smatch match;
   
    std::ifstream infile("./wave_packet_20_20_20_0.1.dat");
    std::string line;
   
    while (std::getline(infile, line)) {
      if (std::regex_search(line, match, rgx)) {
        real=std::atof(match[1].str().c_str());
        imag=std::atof(match[2].str().c_str());
        complex cin(real,imag);
        vecin.push_back(cin);
      }
    }
  for (int i=0;i<norb;++i){
    psi0.set(i,vecin[i]);
    }
  }
  else{
    // Load previous psi values if not first timestep
    PetscViewer    viewer_load;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"./psi_temp/psi-bin.dat",FILE_MODE_READ,&viewer_load);
    VecCreate(PETSC_COMM_WORLD,&psi0.V);
    VecLoad(psi0.V,viewer_load);
    PetscViewerDestroy(&viewer_load);
  }

  get_diffusion(dataa, coord, psi0);  
  save_charge(dataa, coord, dataa.qs, psi0);

  std::vector<Eigen::Vector3d> forces0(dataa.nn, Eigen::Vector3d(3));
  forces0 = force_model(dataa, coord, psi0, iorbs);  
  savetxt("forces0.dat", forces0);
 
  Psidynamics<ETRS_static_stepper, Simple_monitor> pd(psi0, dt, nt, stepper, monitor );
  pd.run();

  // Save current psi values
  PetscViewer    viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"./psi_temp/psi-bin.dat",FILE_MODE_WRITE,&viewer);
  VecView(pd.psi.V,viewer);
  PetscViewerDestroy(&viewer);

  auto forces1 = force_model(dataa, coord, pd.psi,iorbs);  

  // To Convert from Newtons to kCal/Mol/A 
  double N2eVang = 6.241509074E+8;

  for (int j=0;j<nlocal;++j){
      // "id" starts at 1 
      Eigen::Vector3d fpsi { forces1[(id[j]-1)][0], forces1[(id[j]-1)][1], forces1[(id[j]-1)][2]};
      for (int k=0;k<3;++k){ 
          f[j][k]=(fpsi[k]*N2eVang);
      }
  }

  delete [] xx;
  
}
