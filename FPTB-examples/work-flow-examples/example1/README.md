Shows how to get a tight binding band structure from DFT + Wannierization

1) run scf and nscf calculations using pw.x
2) run pw3wannier90.x
3) run wannier90.x
4) delete the first few lines of *_hr.dat, leaving the part of the file that is in a 7-column format, save this as hopping0.dat
5) python filter_hops.py to remove hermitian conjugate matrix elements from hopping0.dat
6) julia bands.jl to compute the tight binding bands
7) python plot.py to plot the tight binding bands
