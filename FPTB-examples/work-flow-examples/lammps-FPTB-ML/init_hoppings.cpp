#include "perovskite_tb.cpp"
#include "blot/lookuptable.hpp"

struct Hoplut{

  std::vector<std::vector<int>> pbp;
  std::vector<std::vector<int>> pbs;
  std::vector<std::vector<int>> brp_par;
  std::vector<std::vector<int>> brp_perp;
  std::vector<std::vector<int>> spsig;
  std::vector<std::vector<int>> ppsig;
  std::vector<std::vector<int>> pppi;
  std::vector<std::vector<int>> spn_pb_to_br;
  std::vector<std::vector<int>> ppn_pb_to_br;
  std::vector<std::vector<int>> ppn_br_to_pb;
  std::vector<std::vector<int>> socpb;
  std::vector<std::vector<int>> socbr_par;
  std::vector<std::vector<int>> socbr_perp;

  Hoplut(){};

  void setup(const PerovData & data){



    //
    //set hoppings
    //
    
    int hopcount = 0;

    //onsite terms
    for (auto iorb : data.orblut.lookup_idx(Orb_typ("Pb:p"))){
      std::vector<int> h = {0, 0, 0,iorb.value, iorb.value};
      pbp.push_back(h);
      hopcount = hopcount+1;
    }

    for (auto iorb : data.orblut.lookup_idx(Orb_typ("Pb:s"))){
      std::vector<int> h = {0, 0, 0, iorb.value, iorb.value};
      pbs.push_back(h);
      hopcount = hopcount+1;
    }

    for (auto iorb : data.orblut.lookup_idx(Orb_typ("Br:p"))){
      auto ddir = data.orblut.lookup_col<Ddir>(iorb);
      auto pdir = data.orblut.lookup_col<Pdir>(iorb);
      if (ddir==pdir){
        std::vector<int> h = { 0, 0, 0,iorb.value, iorb.value};
        brp_par.push_back(h);
        hopcount = hopcount+1;
      }
      else{
        std::vector<int> h = { 0, 0, 0,iorb.value, iorb.value};
        brp_perp.push_back(h);
        hopcount = hopcount+1;
      }
    }

    //nearest neighbor terms
    for (auto ix=0; ix<data.nx; ++ix){
      for (auto iy=0; iy<data.ny; ++iy){
        for (auto iz=0; iz<data.nz; ++iz){

          auto const atom1 = *data.atomlut.lookup_idx(Atomehr_typ("Pb"),Atomehr_ix(ix+1),Atomehr_iy(iy+1),Atomehr_iz(iz+1)).begin();

          for(auto ii : data.shelllut.lookup_idx(Atomehr(atom1),Shell("pbbr_shell_1"))){
            auto atom2 = data.shelllut.lookup_col<Atomehr_nbr>(Shell_idx(ii));
            auto dx = data.shelllut.lookup_col<Dx>(Shell_idx(ii));
            auto dy = data.shelllut.lookup_col<Dy>(Shell_idx(ii));
            auto dz = data.shelllut.lookup_col<Dz>(Shell_idx(ii));
            

            for( auto orb1 : data.orblut.lookup_idx(Atomehr(atom1))){
              for( auto orb2 : data.orblut.lookup_idx(Atomehr(atom2.value))){

                auto ddir2 = data.orblut.lookup_col<Ddir>(orb2);

                auto pdir1 = data.orblut.lookup_col<Pdir>(orb1);
                auto pdir2 = data.orblut.lookup_col<Pdir>(orb2);

                auto spin1 = data.orblut.lookup_col<Spin>(orb1);
                auto spin2 = data.orblut.lookup_col<Spin>(orb2);
                if (data.orblut.lookup_col<Orb_typ>(orb1) == "Pb:s"){
                  if (ddir2==pdir2 && spin1==spin2){
                    std::vector<int> h = {dx.value, dy.value, dz.value, orb1.value, orb2.value};
                    spsig.push_back(h);
                    hopcount = hopcount+1;
                  }
                  else if(spin1==spin2){
                    std::vector<int> h = {dx.value, dy.value, dz.value, orb1.value, orb2.value};
                    //ppsig.push_back(h);
                    spn_pb_to_br.push_back(h);
                    hopcount = hopcount+1;
                  }
                }
                else if (data.orblut.lookup_col<Orb_typ>(orb1) == "Pb:p"){
                  if (ddir2==pdir1 && spin1==spin2){
                    if (ddir2==pdir2){
                      std::vector<int> h = {dx.value, dy.value, dz.value, orb1.value, orb2.value};
                      ppsig.push_back(h);
                      hopcount = hopcount+1;
                    }
                    else{
                      std::vector<int> h = {dx.value, dy.value, dz.value, orb1.value, orb2.value};
                      ppn_pb_to_br.push_back(h);
                      hopcount = hopcount+1;
                    }
                  }
                  else if (ddir2!=pdir1 && spin1==spin2){
                    if (pdir1==pdir2){
                      std::vector<int> h = {dx.value, dy.value, dz.value, orb1.value, orb2.value};
                      pppi.push_back(h);
                      hopcount = hopcount+1;
                    }
                    else if (ddir2==pdir2){
                      std::vector<int> h = {-dx.value, -dy.value, -dz.value, orb2.value, orb1.value};
                      ppn_br_to_pb.push_back(h);
                      hopcount = hopcount+1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  //spin orbit terms
  for (auto ix=0; ix<data.nx; ++ix){
    for (auto iy=0; iy<data.ny; ++iy){
      for (auto iz=0; iz<data.nz; ++iz){

        for (auto atom1 : data.atomlut.lookup_idx(Atomehr_typ("Pb"), Atomehr_ix(ix+1), Atomehr_iy(iy+1), Atomehr_iz(iz+1)) ){

          int orb1 = data.pbp_idx(Atomehr(atom1), Spin("up"), Pdir("x") );
          int orb2 = data.pbp_idx(Atomehr(atom1), Spin("dn"), Pdir("x") );
          int orb3 = data.pbp_idx(Atomehr(atom1), Spin("up"), Pdir("y") );
          int orb4 = data.pbp_idx(Atomehr(atom1), Spin("dn"), Pdir("y") );
          int orb5 = data.pbp_idx(Atomehr(atom1), Spin("up"), Pdir("z") );
          int orb6 = data.pbp_idx(Atomehr(atom1), Spin("dn"), Pdir("z") );

          std::vector<int> h1 = { 0, 0, 0, orb1, orb3}; socpb.push_back(h1); hopcount++;
          std::vector<int> h2 = { 0, 0, 0, orb2, orb4}; socpb.push_back(h2); hopcount++;
          std::vector<int> h3 = { 0, 0, 0, orb2, orb5}; socpb.push_back(h3); hopcount++;
          std::vector<int> h4 = { 0, 0, 0, orb1, orb6}; socpb.push_back(h4); hopcount++;
          std::vector<int> h5 = { 0, 0, 0, orb3, orb6}; socpb.push_back(h5); hopcount++;
          std::vector<int> h6 = { 0, 0, 0, orb4, orb5}; socpb.push_back(h6); hopcount++;
        }

        for (auto atom1 : data.atomlut.lookup_idx(Atomehr_typ("Br1"), Atomehr_ix(ix+1), Atomehr_iy(iy+1), Atomehr_iz(iz+1)) ){

          int orb1 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("x") );
          int orb2 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("x") );
          int orb3 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("y") );
          int orb4 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("y") );
          int orb5 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("z") );
          int orb6 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("z") );

          std::vector<int> h1 = {0, 0, 0, orb1, orb3}; socbr_par.push_back(h1); hopcount++;
          std::vector<int> h2 = {0, 0, 0, orb2, orb4}; socbr_par.push_back(h2); hopcount++;
          std::vector<int> h3 = {0, 0, 0, orb2, orb5}; socbr_par.push_back(h3); hopcount++;
          std::vector<int> h4 = {0, 0, 0, orb1, orb6}; socbr_par.push_back(h4); hopcount++;
          std::vector<int> h5 = {0, 0, 0, orb3, orb6}; socbr_perp.push_back(h5); hopcount++;
          std::vector<int> h6 = {0, 0, 0, orb4, orb5}; socbr_perp.push_back(h6); hopcount++;
        }

        for (auto atom1 : data.atomlut.lookup_idx(Atomehr_typ("Br2"), Atomehr_ix(ix+1), Atomehr_iy(iy+1), Atomehr_iz(iz+1)) ){

          int orb1 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("x") );
          int orb2 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("x") );
          int orb3 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("y") );
          int orb4 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("y") );
          int orb5 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("z") );
          int orb6 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("z") );

          std::vector<int> h1 = {0, 0, 0, orb1, orb3}; socbr_par.push_back(h1); hopcount++;
          std::vector<int> h2 = {0, 0, 0, orb2, orb4}; socbr_par.push_back(h2); hopcount++;
          std::vector<int> h3 = {0, 0, 0, orb2, orb5}; socbr_perp.push_back(h3); hopcount++;
          std::vector<int> h4 = {0, 0, 0, orb1, orb6}; socbr_perp.push_back(h4); hopcount++;
          std::vector<int> h5 = {0, 0, 0, orb3, orb6}; socbr_par.push_back(h5); hopcount++;
          std::vector<int> h6 = {0, 0, 0, orb4, orb5}; socbr_par.push_back(h6); hopcount++;
        }

        for (auto atom1 : data.atomlut.lookup_idx(Atomehr_typ("Br3"), Atomehr_ix(ix+1), Atomehr_iy(iy+1), Atomehr_iz(iz+1)) ){

          int orb1 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("x") );
          int orb2 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("x") );
          int orb3 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("y") );
          int orb4 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("y") );
          int orb5 = data.brp_idx(Atomehr(atom1), Spin("up"), Pdir("z") );
          int orb6 = data.brp_idx(Atomehr(atom1), Spin("dn"), Pdir("z") );

          std::vector<int> h1 = {0, 0, 0, orb1, orb3}; socbr_perp.push_back(h1); hopcount++;
          std::vector<int> h2 = {0, 0, 0, orb2, orb4}; socbr_perp.push_back(h2); hopcount++;
          std::vector<int> h3 = {0, 0, 0, orb2, orb5}; socbr_par.push_back(h3); hopcount++;
          std::vector<int> h4 = {0, 0, 0, orb1, orb6}; socbr_par.push_back(h4); hopcount++;
          std::vector<int> h5 = {0, 0, 0, orb3, orb6}; socbr_par.push_back(h5); hopcount++;
          std::vector<int> h6 = {0, 0, 0, orb4, orb5}; socbr_par.push_back(h6); hopcount++;
         }
        }
      }
    }
  };
};
