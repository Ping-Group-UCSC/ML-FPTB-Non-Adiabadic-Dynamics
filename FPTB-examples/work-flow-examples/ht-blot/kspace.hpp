#ifndef KSPACE
#define KSPACE

#include "blot/blotmath.hpp"

//kpts[i] = ith kpoint
//n = number of samples along each segment, endpoint exclusive
vector<vector<double> > kpath(vector<vector<double> > kpts, int n){
  vector<vector<double> > ans;
  auto m = kpts.size();
  auto p = kpts[0];
  for (int i=0; i<m-1; ++i){
    auto dp = (kpts[i+1] - kpts[i]) / (1.0*n);
    for (int j=0; j<n; ++j){
      ans.push_back(p);
      p = p+dp;
    }
  }
  return ans;
}
      


#endif
