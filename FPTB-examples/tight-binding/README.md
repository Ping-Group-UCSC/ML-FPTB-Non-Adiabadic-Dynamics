General purpose tight binding code
===

bands.jl calculates the band structure of a tight-binding model, given a set of hopping parameters.

How to use
===

The code reads the input file "hopping.dat", which is a list of hopping parameters. This file contains one hopping parameter on each line, with the format

dx dy dz orb1 orb2 Re(h) Im(h)

where dx dy dz is the integer displacement between the two orbitals in units of supercell lattice vectors.

For example: for hopping within the same supercell, dx dy dz = 0 0 0. For hopping from one supercell to the next supercell in the y direction, dx dy dz = 0 1 0. 

orb1 and orb2 are indexes of orbitals 1 and 2.
Re(h) and Im(h) are the real and imaginary parts of the hopping parameter.

The total number of orbitals, nh, is required to be specified in the code. For halide perovskites with a nx*ny*nz supercell, nh = 26*nx*ny*nz. 

By default, the code prints the band eigenvalues to the "work" folder. 