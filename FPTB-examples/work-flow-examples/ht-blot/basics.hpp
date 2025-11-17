#ifndef BASICS
#define BASICS

#include <vector>
#include <memory>
#include <algorithm>


//indexing

//index sequence 
std::vector<int> idxseq(int n1){
  std::vector<int> ans;
  for (auto ii=0;ii<n1;++ii){
    ans.push_back(ii);
  }
  return ans;
}

//endpoint exclusive 
std::vector<int> idxseq(int start, int stop1 ){
  std::vector<int> ans;
  for (auto ii=start;ii<stop1;++ii){
    ans.push_back(ii);
  }
  return ans;
}

//endpoint exclusive 
std::vector<int> idxseq(int start, int stop1, int stride){
  std::vector<int> ans;
  for (auto ii=start;ii<stop1; ii+=stride){
    ans.push_back(ii);
  }
  return ans;
}

struct allseq{};



//sorting

template<class T>
std::vector<int> sortperm(const std::vector<T> & v){
  std::vector<int> idxs;
  for (int i=0; i<v.size(); i++){
    idxs.push_back(i);
  }

  auto comp = [&](const int & a, const int & b){
    return v[a] < v[b];
  };
  sort(idxs.begin(), idxs.end(), comp);

  return idxs;
}

#endif
