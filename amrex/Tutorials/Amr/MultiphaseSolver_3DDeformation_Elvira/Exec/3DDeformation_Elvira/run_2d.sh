rm -rf main2d.gnu.MPI.ex
rm -rf plt*
make DIM=2
mpirun -np 4 ./main2d.gnu.MPI.ex inputs
