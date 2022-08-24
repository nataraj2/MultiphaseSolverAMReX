# Mutiphase flow solver in AMReX

# Installation and compilation 
```
1. git clone https://github.com/nataraj2/MultiphaseSolverAMReX.git
2. cd MultiphaseSolverAMReX
3. cd interfacereconstructionlibrary
4. mkdir build
5. cd build
6. cmake -D EIGEN_INCLUDE_DIR=../eigen -D CMAKE_INSTALL_PREFIX=../install ..
7. make
8. make install
9. cd ../../amrex
10. vi Tools/GNUMake/Make.local
```
In Make.local, change the path to IRLDIR, and the compilers CXX, CC, FC, F90. 
Make sure to use MPI compilers.

# Running
The following is to run the case of Zalesak's disk.
```
1. cd Tutorials/Amr/MultiphaseSolver_Zalesak_Elvira/Exec/Zalesak_Elvira
2. make -j
3. sh run_3d.sh
```
## Input file description
1. `geometry.prob_lo` and `geometry.prob_hi` define the domain
2. `amr.max_level` is the number of refinement levels above the base level. It is set to 2. 
Hence there will be a total of 3 levels (including base level).
3. `amr.regrid_int` is the frequency of regrdding. Has to be 1. Regridding has to be 
done at every step since we do not allow for the interface to cut across levels. The interface 
should always be enclosed by the finest level.

# Visualization
VisIt 2.13 is the minimum version needed for visualization. The AMReX files cannot 
be read by the earlier versions of VisIt. In the running directory, do
```
sh run_output.sh
```
will create a text file `movie.visit` (which contains a list of headers of the files which has all 
the flow variables). It will create a directory `ensight-3D`, in which all plic data is 
written. A file `plic.visit` is also written inside `ensight-3D`. Both `movie.visit` and `plic.visit`
 can be loaded simultaneously into VisIt and visualized.


