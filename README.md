# Mutiphase flow solver in AMReX

# Installation and compilation 
1. ```git clone https://github.com/nataraj2/MultiphaseSolverAMReX.git``
2. ```cd MultiphaseSolverAMReX/```
3. ```cd interfacereconstructionlibrary```
mkdir build
cd build
cmake -D EIGEN_INCLUDE_DIR=../eigen -D CMAKE_INSTALL_PREFIX=../install ..
make
make install
cd ../../amrex
vi Tools/GNUMake/Make.local - Change path to IRLDIR, and the compilers CXX, CC, FC, F90

