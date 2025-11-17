Prerequisites

tested with gcc 9.3.0, PETSc/SLEPc 3.16.0, Eigen 3.4.0  

Requires std=c++2a  

Notes on Eigen:  
Always pass Eigen objects by reference
Don't use auto with Eigen objects
We are compiling with std=c++2a, which avoids alignment issues  
https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html  

Notes on Petsc wrappers:
Recommended to pass VecWrap and MatWrap by reference.
