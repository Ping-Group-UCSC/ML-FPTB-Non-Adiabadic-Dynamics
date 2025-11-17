#ifndef BLOTIO
#define BLOTIO

#include <slepceps.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

template<class T>
vector<vector<T> > loadtxt(const string filename){
  vector<vector<T> > ans;
  fstream ifs;

  ifs.open(filename);
  if (ifs.fail()) cout << "read file error " << filename << endl;

  while (true){
    string line;
    T buf;
    getline(ifs, line);

    stringstream ss(line, ios_base::out|ios_base::in|ios_base::binary);

    if(!ifs) break;

    if (line[0] == '#' || line.empty()) continue;

    vector<T> row;

    while (ss >> buf ) row.push_back(buf);

    ans.push_back(row);
  }

  ifs.close();

  return ans;
}

template<class T>
vector<T> loadtxtvec(const string filename){
  vector<T> ans;
  fstream ifs;

  ifs.open(filename);
  if (ifs.fail()) cout << "read file error " << filename << endl;

  while (true){
    string line;
    T buf;
    getline(ifs, line);

    stringstream ss(line, ios_base::out|ios_base::in|ios_base::binary);

    if(!ifs) break;

    if (line[0] == '#' || line.empty()) continue;

    ss >> buf ;

    ans.push_back(buf);
  }

  ifs.close();

  return ans;
}



template<class T>
void savetxt(const string filename,vector<vector<T> > m){
  ofstream ofs;
  ofs.open(filename);
  for(auto r : m){
    for (auto x : r){
      ofs << x << ' ' ;
    }
    ofs << endl;
  }
  ofs.close();
}

template<class T>
void savetxt(const string filename,vector<T>  m){
  ofstream ofs;
  ofs.open(filename);
  for (auto x : m){
    ofs << x << endl ;
  }
  ofs.close();
}

template<class T>
void savetxt(const string filename, T * m, int size){
  ofstream ofs;
  ofs.open(filename);
  for (int i=0; i<size; ++i){
    ofs << m[i] << endl ;
  }
  ofs.close();
}
 
void savetxt(const string filename, Vec v){//Petsc Vec 
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename.c_str(), &viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(v,viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);
}


#endif
