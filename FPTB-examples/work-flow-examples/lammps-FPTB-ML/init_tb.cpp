#include "perovskite_tb.cpp"
#include "blot/lookuptable.hpp"

PerovData init(){


  int nx = 20;
  int ny = 20;
  int nz = 20;
  int nn = 5*nx*ny*nz;
  double lx = 117.91544768262246;
  double ly = 117.0687848976395 ;
  double lz = 117.17259762488408;
  double lx0 = lx/nx;
  double ly0 = ly/ny;
  double lz0 = lz/nz;

  Eigen::MatrixXd cell (3,3);
  cell(0,0) = lx;
  cell(1,1) = ly;
  cell(2,2) = lz;

  int mshells_max = 1;

  int vcut_x = 2;
  int vcut_y = 2;
  int vcut_z = 2;

  //
  //set lookuptables
  //

  Atomehrlut atomlut;
  auto atomindexdat = loadtxt<string> ("data_to_load/atom_index.dat");
  for (auto line : atomindexdat){
    int ia = std::stoi(line[0]);
    string atom_typ =  line[1];
    int ix = std::stoi(line[2]);
    int iy = std::stoi(line[3]);
    int iz = std::stoi(line[4]);

    atomlut.insert(Atomehr(ia-1),Atomehr_typ(atom_typ),Atomehr_ix(ix),Atomehr_iy(iy),Atomehr_iz(iz)); //atoms 0-indexed in the lookuptable
  }

  Orblut orblut;
  auto orbindexdat = loadtxt<string> ("data_to_load/orb_index.dat");
  int norb = orbindexdat.size();
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

    const auto atom_idx = *atomlut.lookup_idx(Atomehr_typ(atom_typ),Atomehr_ix(ix),Atomehr_iy(iy),Atomehr_iz(iz)).begin();
    orblut.insert(Orb(iorb),Atomehr(atom_idx));
  }


  //
  //set shells
  //

  //auto pbpb_shells = loadtxt<int> ("data_to_load/shells/pbpb_shells.dat");
  auto pbbr_shells = loadtxt<int> ("data_to_load/shells/pbi_shells.dat");
  //auto pbcs_shells = loadtxt<int> ("data_to_load/shells/pbma_shells.dat");
  //auto br1pb_shells = loadtxt<int> ("data_to_load/shells/i1pb_shells.dat");
  //auto br2pb_shells = loadtxt<int> ("data_to_load/shells/i2pb_shells.dat");
  //auto br3pb_shells = loadtxt<int>("data_to_load/shells/i3pb_shells.dat");
  //auto br1cs_shells = loadtxt<int>("data_to_load/shells/i1ma_shells.dat");
  //auto br2cs_shells = loadtxt<int>("data_to_load/shells/i2ma_shells.dat");
  //auto br3cs_shells = loadtxt<int>("data_to_load/shells/i3ma_shells.dat");
  //auto br1br_shells = loadtxt<int>("data_to_load/shells/i1i_shells.dat");
  //auto br2br_shells = loadtxt<int>("data_to_load/shells/i2i_shells.dat");
  //auto br3br_shells = loadtxt<int>("data_to_load/shells/i3i_shells.dat");

  //std::vector<
  //  std::vector<std::vector<int>>> brpb_shells {br1pb_shells, br2pb_shells,br3pb_shells};
  //std::vector<
  //  std::vector<std::vector<int>>> brcs_shells {br1cs_shells,br2cs_shells,br3cs_shells};
  //std::vector<
  //  std::vector<std::vector<int>>> brbr_shells {br1br_shells,br2br_shells,br3br_shells};

  //auto n_pbpb_shells = loadtxtvec<int> ("data_to_load/shells/n_pbpb_shells.dat");
  auto n_pbbr_shells = loadtxtvec<int> ("data_to_load/shells/n_pbi_shells.dat");
  //auto n_pbcs_shells = loadtxtvec<int> ("data_to_load/shells/n_pbma_shells.dat");
  //auto n_br1pb_shells = loadtxtvec<int> ("data_to_load/shells/n_i1pb_shells.dat");
  //auto n_br2pb_shells = loadtxtvec<int> ("data_to_load/shells/n_i2pb_shells.dat");
  //auto n_br3pb_shells = loadtxtvec<int>("data_to_load/shells/n_i3pb_shells.dat");
  //auto n_br1cs_shells = loadtxtvec<int>("data_to_load/shells/n_i1ma_shells.dat");
  //auto n_br2cs_shells = loadtxtvec<int>("data_to_load/shells/n_i2ma_shells.dat");
  //auto n_br3cs_shells = loadtxtvec<int>("data_to_load/shells/n_i3ma_shells.dat");
  //auto n_br1br_shells = loadtxtvec<int>("data_to_load/shells/n_i1i_shells.dat");
  //auto n_br2br_shells = loadtxtvec<int>("data_to_load/shells/n_i2i_shells.dat");
  //auto n_br3br_shells = loadtxtvec<int>("data_to_load/shells/n_i3i_shells.dat");

  //std::vector<std::vector<int>> n_brpb_shells {n_br1pb_shells, n_br2pb_shells,n_br3pb_shells};
  //std::vector<std::vector<int>> n_brcs_shells {n_br1cs_shells,n_br2cs_shells,n_br3cs_shells};
  //std::vector<std::vector<int>> n_brbr_shells {n_br1br_shells,n_br2br_shells,n_br3br_shells};


  Shelllut shelllut;
  int count = 0;
  for (auto ix=0; ix<nx; ++ix){
    for (auto iy=0; iy<ny; ++iy){
      for (auto iz=0; iz<nz; ++iz){

        auto const atom1 = *atomlut.lookup_idx(Atomehr_typ("Pb"),Atomehr_ix(ix+1),Atomehr_iy(iy+1),Atomehr_iz(iz+1)).begin();

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
            auto  const atom2 = *atomlut.lookup_idx(Atomehr_typ(brtype),Atomehr_ix(dvx.rem+1),Atomehr_iy(dvy.rem+1),Atomehr_iz(dvz.rem+1)).begin();
            
            string shellname = "pbbr_shell_" + std::to_string(ishell+1);

            shelllut.insert(Shell_idx(count),Shell(shellname),Atomehr(atom1),Atomehr_nbr(atom2.value),Dx(dvx.quot),Dy(dvy.quot),Dz(dvz.quot));

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

  //Vcutlut vcutlut;
  Vcutlut vcutlut(nn, std::vector<std::vector<int>>(624,std::vector<int>(5)));
  if (FILE *file = fopen("vcutlut_vector.dat", "r")) {
    fclose(file);
    cout << "Vcut Loaded" << endl;
    //vcutlut = loadtxtveclut<int>("vcutlut_vector.dat");
    std::vector<std::vector<int>> vcutlut_t;
    vcutlut_t = loadtxt<int>("vcutlut_vector.dat");
    for (int i =0; i<nn;++i){
      std::vector<std::vector<int>> nbr_list(624);

      for (int j=0;j<624;++j){
        nbr_list[j] = vcutlut_t[(i*624)+j];
      }
      vcutlut[i] = nbr_list;
    } 
  }
  else{
    for (int ia1=0; ia1<nn; ++ia1){
      Atomehr atom1(ia1);

      std::vector<std::vector<int>> nbr_list;

      auto ix = atomlut.lookup_col<Atomehr_ix>(atom1);
      auto iy = atomlut.lookup_col<Atomehr_iy>(atom1);
      auto iz = atomlut.lookup_col<Atomehr_iz>(atom1);

      for (auto dx=-vcut_x; dx<vcut_x+1; ++dx){
        for (auto dy=-vcut_y; dy<vcut_y+1; ++dy){
          for (auto dz=-vcut_z; dz<vcut_z+1; ++dz){

            auto dvx = knuthdiv( ix.value-1 + dx , nx);
            auto dvy = knuthdiv( iy.value-1 + dy , ny);
            auto dvz = knuthdiv( iz.value-1 + dz , nz);

            for( auto atom2 : atomlut.lookup_idx(Atomehr_ix(dvx.rem+1), Atomehr_iy(dvy.rem+1), Atomehr_iz(dvz.rem+1))){

              if (atom1==atom2 && dx==0 && dy==0 && dz==0) continue;
              std::vector<int> nbr {atom1.value, atom2.value, dvx.quot, dvy.quot, dvz.quot};
              nbr_list.push_back(nbr);
            }
          }
        }
      }
      vcutlut[ia1] = nbr_list;
    }
    savetxtveclut("vcutlut_vector.dat",vcutlut);
  }

  PerovData ans(atomlut,orblut,shelllut, vcutlut, nx,ny,nz,nn,norb, cell);
  return ans;
}

