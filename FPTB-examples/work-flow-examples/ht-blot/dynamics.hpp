#ifndef DYNAMICS
#define DYNAMICS

#include <slepceps.h>
#include <petscts.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>
#include <memory>
#include "blot/blotio.hpp"
#include "blot/blotmath.hpp"
#include "blot/basics.hpp"
#include "blot/units.hpp"



//class for running dynamics specified by a time-dependent hamiltonian.
//To use, specify the following:
// a functor Stepper_t // this performs one time step, changing psi into psi1
// a functor Monitor_t // for printing the trajectory
template<class Stepper_t, class Monitor_t>
struct Psidynamics{
  VecWrap& psi;
  double dt;
  int nt; //total number of time steps
  int iit; //current time step

  Stepper_t& stepper;
  Monitor_t& monitor;

  void run();

  Psidynamics(VecWrap& psi0, double dt0, int nt0, Stepper_t& stepper0, Monitor_t& monitor0);
};


template<class Stepper_t, class Monitor_t>
Psidynamics<Stepper_t, Monitor_t>::Psidynamics(VecWrap& psi0, double dt0, int nt0, Stepper_t& stepper0, Monitor_t& monitor0) :
psi(psi0), dt(dt0), nt(nt0), stepper(stepper0), monitor(monitor0) { 
  iit=0;
  //stepper.init(psi, dt, nt, iit);
}

template<class Stepper_t, class Monitor_t>
void Psidynamics<Stepper_t, Monitor_t>::run(){
  VecWrap psi1(psi.n);

  while (iit<nt){
    stepper(psi, psi1, dt, nt, iit);
    VecCopy(psi1.V, psi.V);
    monitor(psi, dt, nt, iit);
    ++iit;
  }
}

//ETRS timestepper

struct ETRS_static_stepper{
  MatWrap & h;
  ETRS_static_stepper(MatWrap & h0) : h(h0) {};
  void operator() (VecWrap & psi, VecWrap & psi1, double dt, int nt, int iit);
};

// d/dt (psi) = (i/hbar) H psi
// units of H = eV
// units of t = ps
void ETRS_static_stepper::operator() (VecWrap& psi, VecWrap& psi1, double dt, int nt, int iit){
  PetscErrorCode ierr;

  VecWrap psi_temp(psi.n);
  VecWrap psi_temp1(psi.n);
  VecCopy(psi.V, psi_temp.V);

  //0th order
  VecCopy(psi_temp.V, psi1.V);

  for (auto jj=1; jj<4; ++jj){
    MatMult( h.A, psi_temp.V, psi_temp1.V );
    VecScale(psi_temp1.V,  dt*(1.0/jj)* PETSC_i * elemSI/hbarSI/1e12);
    VecCopy(psi_temp1.V, psi_temp.V);
    VecAXPY(psi1.V, 1.0, psi_temp.V);
  }
}

struct Dummy_monitor{
  Dummy_monitor(int) {};
  void operator() (VecWrap& psi, double dt, int nt, int iit) {};
};

struct Simple_monitor{
  int stepsize; //stride for saving the solution
  string filename;
  std::vector<int> orbitals; //which orbitals to save, 0-indexed
  ofstream ofs;

  Simple_monitor(int stepsize1, string filename1, std::vector<int> orbitals1);
  ~Simple_monitor();

  void operator() (VecWrap& psi, double dt, int nt, int iit);
};

Simple_monitor::Simple_monitor(int stepsize1, string filename1, std::vector<int> orbitals1):
stepsize(stepsize1) , filename(filename1), orbitals(orbitals1){
  ofs.open(filename);
}

Simple_monitor::~Simple_monitor(){
  ofs.close();
}

void Simple_monitor::operator() (VecWrap& psi, double dt, int nt, int iit){
  
  if (iit % stepsize == 0){
    //VecGetValues(psi.V, orbitals.n, orbitals.array(), to_save.get());
    
    PetscScalar *arr;
    VecGetArray(psi.V, &arr);
    for(auto ii=0; ii < orbitals.size(); ++ii){
      ofs << arr[orbitals[ii]] << " " ;
    }
    ofs << endl;
    VecRestoreArray(psi.V, &arr);

  }
  
}
 


//
// Use Petsc's timestepping library
// 
//

extern "C" typedef PetscErrorCode RHS_func(TS ts, PetscReal t, Vec u, Vec F, void * h_c);
extern "C" typedef int Monitor_func(TS ts,PetscInt steps,PetscReal time,Vec u,void * m);

struct Psidynamics_PetscTS{
  VecWrap psi;
  double dt;
  int nt;

  RHS_func * h_f;
  void * h_c;
  
  Monitor_func * m_f;
  void * m_c;

  TS ts;

  PetscErrorCode run();
  Psidynamics_PetscTS(VecWrap psi0, RHS_func * h_f0, void * h_c0, Monitor_func * m_f0, void * m_c0, double dt0, int nt0);
};

Psidynamics_PetscTS::Psidynamics_PetscTS(VecWrap psi0, RHS_func * h_f0, void * h_c0, Monitor_func * m_f0, void * m_c0, double dt0, int nt0) :
psi(psi0), h_f(h_f0), h_c(h_c0), m_f(m_f0), m_c(m_c0), dt(dt0), nt(nt0) {

  TSCreate(PETSC_COMM_WORLD,&ts);
  TSSetProblemType(ts,TS_LINEAR);
  TSSetSolution(ts,psi.V);
  TSSetType(ts,TSRK);
  //TSSetType(ts,TSSSP);

  TSSetRHSFunction(ts, NULL, h_f, h_c);

  TSSetTimeStep(ts,dt);
  TSSetMaxSteps(ts,nt);
  TSSetMaxTime(ts,nt*dt);

  TSMonitorSet(ts, m_f, m_c, NULL);

}

// d/dt (psi) = (i/hbar) H psi
// units of H = eV
// units of t = ps
// MatWrap * h_c is the hamiltonian matrix
extern "C" PetscErrorCode Static_Ham_func(TS ts, PetscReal t, Vec u, Vec F, void * h_c){
  PetscErrorCode ierr;
  auto h_c1 = static_cast<MatWrap*>(h_c);
  
  //ierr = VecCopy(u,F);
  //std::cout << "test" << endl;

  ierr = MatMult( h_c1->A, u, F);
  ierr = VecScale(F, PETSC_i * elemSI/hbarSI/1e12);

  return ierr;
}

PetscErrorCode Psidynamics_PetscTS::run(){
  PetscErrorCode ierr;
  ierr = TSSolve(ts, psi.V);
  return ierr;
}

// functions for monitoring progress of trajectory

struct Simple_monitor_PetscTS{
  int stepsize; //stride for saving the solution
  const string filename;
  std::vector<int> orbitals; //which orbitals to save
  ofstream ofs;
  Simple_monitor_PetscTS(int stepsize1, const string filename1, std::vector<int> orbitals1);
  ~Simple_monitor_PetscTS();
};
 
Simple_monitor_PetscTS::Simple_monitor_PetscTS(int stepsize1, const string filename1, std::vector<int> orbitals1):
stepsize(stepsize1) , filename(filename1), orbitals(orbitals1){
  ofs.open(filename);
}

Simple_monitor_PetscTS::~Simple_monitor_PetscTS(){
  ofs.close();
}

extern "C" int Simple_monitor_PetscTS_func (TS ts,PetscInt steps,PetscReal time,Vec u,void * m_c){
  
  Simple_monitor_PetscTS* m_c1 = static_cast<Simple_monitor_PetscTS*>(m_c);
  if (steps % m_c1->stepsize == 0){

    std::shared_ptr<PetscScalar[]> to_save( new PetscScalar[m_c1->orbitals.size()]);
    VecGetValues(u, m_c1->orbitals.size(), &(m_c1->orbitals[0]), to_save.get());
    for(auto ii=0; ii < m_c1->orbitals.size(); ++ii){
      m_c1->ofs << to_save.get()[ii] << " " ;
    }
    m_c1->ofs << endl;
  }
 
  return 0;
}

  
#endif
