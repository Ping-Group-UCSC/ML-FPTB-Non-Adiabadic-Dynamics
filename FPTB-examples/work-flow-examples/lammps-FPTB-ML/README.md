Performs time-independent Hamiltonian dynamics with DEEPMD Machine Learning Force Fields 
Uses lammps compiled with DEEPMD support

1. Modify makefile in this directory to match your system

2. Compile `ehrenfest` using makefile located in this directory
  
  make ehrenfest

3. To run use:

        mpirun ./ehrenfest lammps1.in  > log.out


Units:  
[Hamiltonian matrix elements] = eV  
[time] = ps  
[atomic coordinates] = angstrom  
[md_charges] = e  
[forces] = N  
