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
In Make.local, change the path to IRLDIR, and the compilers CXX, CC, FC, F90

# How to run?
```
cd Tutorials/Amr/MultiphaseSolver_Zalesak_Elvira/Exec/Zalesak_Elvira
make -j
sh run_3d.sh
```


