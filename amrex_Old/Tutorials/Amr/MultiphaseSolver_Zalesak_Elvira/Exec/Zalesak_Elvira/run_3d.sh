rm -rf main3d.gnu.MPI.ex
rm -rf chk*
rm -rf plt*
rm -rf ensight-3D
mkdir ensight-3D
make -j DIM=3
#mpirun -np 1 ./main3d.gnu.DEBUG.MPI.ex inputs
mpiexec -np 4 ./main3d.gnu.MPI.ex inputs
#gdb exec-file main3d.gnu.MPI.ex symbol-file inputs
