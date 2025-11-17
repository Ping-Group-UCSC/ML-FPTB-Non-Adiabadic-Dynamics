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

  Hoplut(const PerovData & data);

};

Hoplut::Hoplut(const PerovData & data){

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

        auto const atom1 = *data.atomlut.lookup_idx(Atom_typ("Pb"),Atom_ix(ix+1),Atom_iy(iy+1),Atom_iz(iz+1)).begin();

        for(auto ii : data.shelllut.lookup_idx(Atom(atom1),Shell("pbbr_shell_1"))){
          auto atom2 = data.shelllut.lookup_col<Atom_nbr>(Shell_idx(ii));
          auto dx = data.shelllut.lookup_col<Dx>(Shell_idx(ii));
          auto dy = data.shelllut.lookup_col<Dy>(Shell_idx(ii));
          auto dz = data.shelllut.lookup_col<Dz>(Shell_idx(ii));
          

          for( auto orb1 : data.orblut.lookup_idx(Atom(atom1))){
            for( auto orb2 : data.orblut.lookup_idx(Atom(atom2.value))){

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

        for (auto atom1 : data.atomlut.lookup_idx(Atom_typ("Pb"), Atom_ix(ix+1), Atom_iy(iy+1), Atom_iz(iz+1)) ){

          auto const orb1 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Pb:p"), Pdir("x") ).begin(); 
          auto const orb2 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Pb:p"), Pdir("x") ).begin(); 
          auto const orb3 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Pb:p"), Pdir("y") ).begin(); 
          auto const orb4 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Pb:p"), Pdir("y") ).begin(); 
          auto const orb5 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Pb:p"), Pdir("z") ).begin(); 
          auto const orb6 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Pb:p"), Pdir("z") ).begin(); 

          std::vector<int> h1 = { 0, 0, 0, orb1.value, orb3.value}; socpb.push_back(h1); hopcount++;
          std::vector<int> h2 = { 0, 0, 0, orb2.value, orb4.value}; socpb.push_back(h2); hopcount++;
          std::vector<int> h3 = { 0, 0, 0, orb2.value, orb5.value}; socpb.push_back(h3); hopcount++;
          std::vector<int> h4 = { 0, 0, 0, orb1.value, orb6.value}; socpb.push_back(h4); hopcount++;
          std::vector<int> h5 = { 0, 0, 0, orb3.value, orb6.value}; socpb.push_back(h5); hopcount++;
          std::vector<int> h6 = { 0, 0, 0, orb4.value, orb5.value}; socpb.push_back(h6); hopcount++;
        }

        for (auto atom1 : data.atomlut.lookup_idx(Atom_typ("Br1"), Atom_ix(ix+1), Atom_iy(iy+1), Atom_iz(iz+1)) ){

          auto const orb1 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("x") ).begin(); 
          auto const orb2 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("x") ).begin(); 
          auto const orb3 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("y") ).begin(); 
          auto const orb4 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("y") ).begin(); 
          auto const orb5 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("z") ).begin(); 
          auto const orb6 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("z") ).begin(); 

          std::vector<int> h1 = {0, 0, 0, orb1.value, orb3.value}; socbr_par.push_back(h1); hopcount++;
          std::vector<int> h2 = {0, 0, 0, orb2.value, orb4.value}; socbr_par.push_back(h2); hopcount++;
          std::vector<int> h3 = {0, 0, 0, orb2.value, orb5.value}; socbr_par.push_back(h3); hopcount++;
          std::vector<int> h4 = {0, 0, 0, orb1.value, orb6.value}; socbr_par.push_back(h4); hopcount++;
          std::vector<int> h5 = {0, 0, 0, orb3.value, orb6.value}; socbr_perp.push_back(h5); hopcount++;
          std::vector<int> h6 = {0, 0, 0, orb4.value, orb5.value}; socbr_perp.push_back(h6); hopcount++;
        }

        for (auto atom1 : data.atomlut.lookup_idx(Atom_typ("Br2"), Atom_ix(ix+1), Atom_iy(iy+1), Atom_iz(iz+1)) ){

          auto const orb1 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("x") ).begin(); 
          auto const orb2 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("x") ).begin(); 
          auto const orb3 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("y") ).begin(); 
          auto const orb4 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("y") ).begin(); 
          auto const orb5 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("z") ).begin(); 
          auto const orb6 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("z") ).begin(); 

          std::vector<int> h1 = {0, 0, 0, orb1.value, orb3.value}; socbr_par.push_back(h1); hopcount++;
          std::vector<int> h2 = {0, 0, 0, orb2.value, orb4.value}; socbr_par.push_back(h2); hopcount++;
          std::vector<int> h3 = {0, 0, 0, orb2.value, orb5.value}; socbr_perp.push_back(h3); hopcount++;
          std::vector<int> h4 = {0, 0, 0, orb1.value, orb6.value}; socbr_perp.push_back(h4); hopcount++;
          std::vector<int> h5 = {0, 0, 0, orb3.value, orb6.value}; socbr_par.push_back(h5); hopcount++;
          std::vector<int> h6 = {0, 0, 0, orb4.value, orb5.value}; socbr_par.push_back(h6); hopcount++;
        }

        for (auto atom1 : data.atomlut.lookup_idx(Atom_typ("Br3"), Atom_ix(ix+1), Atom_iy(iy+1), Atom_iz(iz+1)) ){

          auto const orb1 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("x") ).begin(); 
          auto const orb2 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("x") ).begin(); 
          auto const orb3 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("y") ).begin(); 
          auto const orb4 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("y") ).begin(); 
          auto const orb5 = *data.orblut.lookup_idx(Atom(atom1), Spin("up"), Orb_typ("Br:p"), Pdir("z") ).begin(); 
          auto const orb6 = *data.orblut.lookup_idx(Atom(atom1), Spin("dn"), Orb_typ("Br:p"), Pdir("z") ).begin(); 

          std::vector<int> h1 = {0, 0, 0, orb1.value, orb3.value}; socbr_perp.push_back(h1); hopcount++;
          std::vector<int> h2 = {0, 0, 0, orb2.value, orb4.value}; socbr_perp.push_back(h2); hopcount++;
          std::vector<int> h3 = {0, 0, 0, orb2.value, orb5.value}; socbr_par.push_back(h3); hopcount++;
          std::vector<int> h4 = {0, 0, 0, orb1.value, orb6.value}; socbr_par.push_back(h4); hopcount++;
          std::vector<int> h5 = {0, 0, 0, orb3.value, orb6.value}; socbr_par.push_back(h5); hopcount++;
          std::vector<int> h6 = {0, 0, 0, orb4.value, orb5.value}; socbr_par.push_back(h6); hopcount++;
        }
      }
    }
  }

}

