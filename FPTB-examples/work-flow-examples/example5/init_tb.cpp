#include "perovskite_tb.cpp"
#include "blot/lookuptable.hpp"

PerovData init(){
  int nn = 40;
  int nx = 2;
  int ny = 2;
  int nz = 2;
  double lx = 12.00445150246208;
  double ly = 12.00445150246208;
  double lz = 12.00445150246208;
  double lx0 = 12.00445150246208/nx;
  double ly0 = 12.00445150246208/ny;
  double lz0 = 12.00445150246208/nz;

  Eigen::MatrixXd cell (3,3);
  cell(0,0) = lx;
  cell(1,1) = ly;
  cell(2,2) = lz;

  int mshells_max = 1;

  //
  //set lookuptables
  //

  Atomlut atomlut;
  auto atomindexdat = loadtxt<string> ("../example3/atom_index.dat");
  for (auto line : atomindexdat){
    int ia = std::stoi(line[0]);
    string atom_typ =  line[1];
    int ix = std::stoi(line[2]);
    int iy = std::stoi(line[3]);
    int iz = std::stoi(line[4]);

    atomlut.insert(Atom(ia-1),Atom_typ(atom_typ),Atom_ix(ix),Atom_iy(iy),Atom_iz(iz)); //atoms 0-indexed in the lookuptable
  }

  Orblut orblut;
  auto orbindexdat = loadtxt<string> ("../example3/orb_index.dat");
  for (auto line : orbindexdat){
    int iorb = std::stoi(line[0]); //orb index is 1-indexed
    string orb_typ =  line[1];
    int ix = std::stoi(line[2]);
    int iy = std::stoi(line[3]);
    int iz = std::stoi(line[4]);

    string atom_typ;

    if (orb_typ=="Pb:s"){
      string spin = line[5];
      orblut.insert(Orb(iorb),Orb_typ(orb_typ),Spin(spin));
      atom_typ = "Pb";
    }
    else if (orb_typ=="Pb:p"){
      string pdir = line[5];
      string pdir_lower = pdir;
      std::transform(pdir.begin(), pdir.end(), pdir_lower.begin(), (int (*)(int))tolower);
      string spin = line[6];
      orblut.insert(Orb(iorb),Orb_typ(orb_typ),Pdir(pdir_lower),Spin(spin));
      atom_typ = "Pb";
    }
    else if (orb_typ=="Br:p"){
      string ddir = line[5];
      string pdir = line[6];
      string pdir_lower = pdir;
      std::transform(pdir.begin(), pdir.end(), pdir_lower.begin(), (int (*)(int))tolower);
      string spin = line[7];
      orblut.insert(Orb(iorb),Orb_typ(orb_typ),Ddir(ddir),Pdir(pdir_lower),Spin(spin));
      if (ddir=="x") atom_typ = "Br1";
      if (ddir=="y") atom_typ = "Br2";
      if (ddir=="z") atom_typ = "Br3";
    }

    const auto atom_idx = *atomlut.lookup_idx(Atom_typ(atom_typ),Atom_ix(ix),Atom_iy(iy),Atom_iz(iz)).begin();
    orblut.insert(Orb(iorb),Atom(atom_idx));
  }


  //
  //set shells
  //

  //auto pbpb_shells = loadtxt<int> ("../example3/shells/pbpb_shells.dat");
  auto pbbr_shells = loadtxt<int> ("../example3/shells/pbi_shells.dat");
  //auto pbcs_shells = loadtxt<int> ("../example3/shells/pbma_shells.dat");
  //auto br1pb_shells = loadtxt<int> ("../example3/shells/i1pb_shells.dat");
  //auto br2pb_shells = loadtxt<int> ("../example3/shells/i2pb_shells.dat");
  //auto br3pb_shells = loadtxt<int>("../example3/shells/i3pb_shells.dat");
  //auto br1cs_shells = loadtxt<int>("../example3/shells/i1ma_shells.dat");
  //auto br2cs_shells = loadtxt<int>("../example3/shells/i2ma_shells.dat");
  //auto br3cs_shells = loadtxt<int>("../example3/shells/i3ma_shells.dat");
  //auto br1br_shells = loadtxt<int>("../example3/shells/i1i_shells.dat");
  //auto br2br_shells = loadtxt<int>("../example3/shells/i2i_shells.dat");
  //auto br3br_shells = loadtxt<int>("../example3/shells/i3i_shells.dat");

  //std::vector<
  //  std::vector<std::vector<int>>> brpb_shells {br1pb_shells, br2pb_shells,br3pb_shells};
  //std::vector<
  //  std::vector<std::vector<int>>> brcs_shells {br1cs_shells,br2cs_shells,br3cs_shells};
  //std::vector<
  //  std::vector<std::vector<int>>> brbr_shells {br1br_shells,br2br_shells,br3br_shells};

  //auto n_pbpb_shells = loadtxtvec<int> ("../example3/shells/n_pbpb_shells.dat");
  auto n_pbbr_shells = loadtxtvec<int> ("../example3/shells/n_pbi_shells.dat");
  //auto n_pbcs_shells = loadtxtvec<int> ("../example3/shells/n_pbma_shells.dat");
  //auto n_br1pb_shells = loadtxtvec<int> ("../example3/shells/n_i1pb_shells.dat");
  //auto n_br2pb_shells = loadtxtvec<int> ("../example3/shells/n_i2pb_shells.dat");
  //auto n_br3pb_shells = loadtxtvec<int>("../example3/shells/n_i3pb_shells.dat");
  //auto n_br1cs_shells = loadtxtvec<int>("../example3/shells/n_i1ma_shells.dat");
  //auto n_br2cs_shells = loadtxtvec<int>("../example3/shells/n_i2ma_shells.dat");
  //auto n_br3cs_shells = loadtxtvec<int>("../example3/shells/n_i3ma_shells.dat");
  //auto n_br1br_shells = loadtxtvec<int>("../example3/shells/n_i1i_shells.dat");
  //auto n_br2br_shells = loadtxtvec<int>("../example3/shells/n_i2i_shells.dat");
  //auto n_br3br_shells = loadtxtvec<int>("../example3/shells/n_i3i_shells.dat");

  //std::vector<std::vector<int>> n_brpb_shells {n_br1pb_shells, n_br2pb_shells,n_br3pb_shells};
  //std::vector<std::vector<int>> n_brcs_shells {n_br1cs_shells,n_br2cs_shells,n_br3cs_shells};
  //std::vector<std::vector<int>> n_brbr_shells {n_br1br_shells,n_br2br_shells,n_br3br_shells};


  Shelllut shelllut;
  int count = 0;
  for (auto ix=0; ix<nx; ++ix){
    for (auto iy=0; iy<ny; ++iy){
      for (auto iz=0; iz<nz; ++iz){

        auto const atom1 = *atomlut.lookup_idx(Atom_typ("Pb"),Atom_ix(ix+1),Atom_iy(iy+1),Atom_iz(iz+1)).begin();

        /*
        //pbpb shells
        for (auto ishell=0; ishell<mshells_max; ++ishell){
          auto j1shell = n_pbpb_shells[ishell+1]; 
          auto j2shell = n_pbpb_shells[ishell+2];
  
          for (auto jshell=j1shell; jshell<j2shell; ++jshell){
            auto dvx = knuthdiv( ix + pbpb_shells[jshell][0] , nx);
            auto dvy = knuthdiv( iy + pbpb_shells[jshell][1] , ny);
            auto dvz = knuthdiv( iz + pbpb_shells[jshell][2] , nz);
            auto const atom2 = *atomlut.lookup_idx(Atom_typ("Pb"),Atom_ix(dvx.rem+1),Atom_iy(dvy.rem+1),Atom_iz(dvz.rem+1)).begin();
        
  
            string shellname = "pbpb_shell_" + std::to_string(ishell+1);
            shelllut.insert(Shell_idx(count),Shell(shellname),Atom(atom1),Atom_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));
            count = count+1;
          }
        }
        */


        //pbbr shells
        for (auto ishell=0; ishell<mshells_max; ++ishell){
          auto j1shell = n_pbbr_shells[ishell+0]; 
          auto j2shell = n_pbbr_shells[ishell+1];
  
          for (auto jshell=j1shell; jshell<j2shell; ++jshell){
            auto dvx = knuthdiv( ix + pbbr_shells[jshell][0] , nx);
            auto dvy = knuthdiv( iy + pbbr_shells[jshell][1] , ny);
            auto dvz = knuthdiv( iz + pbbr_shells[jshell][2] , nz);
            auto jtyp = pbbr_shells[jshell][3];
            string brtype = "Br" + std::to_string(jtyp);
            auto  const atom2 = *atomlut.lookup_idx(Atom_typ(brtype),Atom_ix(dvx.rem+1),Atom_iy(dvy.rem+1),Atom_iz(dvz.rem+1)).begin();
            
            string shellname = "pbbr_shell_" + std::to_string(ishell+1);

            shelllut.insert(Shell_idx(count),Shell(shellname),Atom(atom1),Atom_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));

            count = count+1;

          }
        }
        
        /*
        //pbcs shells
        for (auto ishell=0; ishell<mshells_max; ++ishell){
          auto j1shell = n_pbcs_shells[ishell+0]; 
          auto j2shell = n_pbcs_shells[ishell+1];
  
          for (auto jshell=j1shell; jshell<j2shell; ++jshell){
            auto dvx = knuthdiv( ix + pbcs_shells[jshell][0] , nx);
            auto dvy = knuthdiv( iy + pbcs_shells[jshell][1] , ny);
            auto dvz = knuthdiv( iz + pbcs_shells[jshell][2] , nz);
            auto const atom2 = *atomlut.lookup_idx(Atom_typ("Cs"),Atom_ix(dvx.rem+1),Atom_iy(dvy.rem+1),Atom_iz(dvz.rem+1)).begin();
  
            string shellname = "pbcs_shell_" + std::to_string(ishell+1);
            shelllut.insert(Shell_idx(count),Shell(shellname),Atom(atom1),Atom_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));
            count = count+1;
          }
        }


        for (auto ityp=1; ityp<4; ++ityp){

          string brtype = "Br" + std::to_string(ityp);
          auto const atom1br = *atomlut.lookup_idx(Atom_typ(brtype),Atom_ix(ix+1),Atom_iy(iy+1),Atom_iz(iz+1)).begin();

          //brpb shells
          for (auto ishell=0; ishell<mshells_max; ++ishell){
            auto j1shell = n_brpb_shells[ityp-1][ishell+0]; 
            auto j2shell = n_brpb_shells[ityp-1][ishell+1];
  
            for (auto jshell=j1shell; jshell<j2shell; ++jshell){
              auto dvx = knuthdiv( ix + brpb_shells[ityp-1][jshell][0] , nx);
              auto dvy = knuthdiv( iy + brpb_shells[ityp-1][jshell][1] , ny);
              auto dvz = knuthdiv( iz + brpb_shells[ityp-1][jshell][2] , nz);
              auto const atom2 = *atomlut.lookup_idx(Atom_typ("Pb"),Atom_ix(dvx.rem+1),Atom_iy(dvy.rem+1),Atom_iz(dvz.rem+1)).begin();
  
              string shellname = "brpb_shell_" + std::to_string(ishell+1);
              shelllut.insert(Shell_idx(count),Shell(shellname),Atom(atom1br),Atom_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));
              count = count+1;
            }
          }


          //brcs shells
          for (auto ishell=0; ishell<mshells_max; ++ishell){
            auto j1shell = n_brcs_shells[ityp-1][ishell+0]; 
            auto j2shell = n_brcs_shells[ityp-1][ishell+1];
  
            for (auto jshell=j1shell; jshell<j2shell; ++jshell){
              auto dvx = knuthdiv( ix + brcs_shells[ityp-1][jshell][0] , nx);
              auto dvy = knuthdiv( iy + brcs_shells[ityp-1][jshell][1] , ny);
              auto dvz = knuthdiv( iz + brcs_shells[ityp-1][jshell][2] , nz);
              auto const atom2 = *atomlut.lookup_idx(Atom_typ("Cs"),Atom_ix(dvx.rem+1),Atom_iy(dvy.rem+1),Atom_iz(dvz.rem+1)).begin();
  
              string shellname = "brcs_shell_" + std::to_string(ishell+1);
              shelllut.insert(Shell_idx(count),Shell(shellname),Atom(atom1br),Atom_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));
              count = count+1;
            }
          }


          //brbr shells
          for (auto ishell=0; ishell<mshells_max; ++ishell){
            auto j1shell = n_brbr_shells[ityp-1][ishell+1]; 
            auto j2shell = n_brbr_shells[ityp-1][ishell+2];
  
            for (auto jshell=j1shell; jshell<j2shell; ++jshell){
              auto dvx = knuthdiv( ix + brbr_shells[ityp-1][jshell][0] , nx);
              auto dvy = knuthdiv( iy + brbr_shells[ityp-1][jshell][1] , ny);
              auto dvz = knuthdiv( iz + brbr_shells[ityp-1][jshell][2] , nz);
              auto jtyp = brbr_shells[ityp-1][jshell][3];
              string brtype = "Br" + std::to_string(jtyp);
              auto const atom2 = *atomlut.lookup_idx(Atom_typ(brtype),Atom_ix(dvx.rem+1),Atom_iy(dvy.rem+1),Atom_iz(dvz.rem+1)).begin();
  
              string shellname = "brbr_shell_" + std::to_string(ishell+1);
              shelllut.insert(Shell_idx(count),Shell(shellname),Atom(atom1br),Atom_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));
              count = count+1;
            }
          }

        }

  */
      }
    }
  }
  

  PerovData ans(atomlut,orblut,shelllut, nx,ny,nz,nn, cell);
  return ans;
}

