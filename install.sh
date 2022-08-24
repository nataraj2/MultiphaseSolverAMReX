cd interfacereconstructionlibrary
mkdir build
cd build
cmake -D EIGEN_INCLUDE_DIR=../eigen -D CMAKE_INSTALL_PREFIX=../install ..
make
make install
