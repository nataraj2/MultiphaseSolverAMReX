# Multiphase flow solver in AMReX
This repository contains the code for an all-Mach multiphase flow solver for liquid-gas 
flows using a geometric volume-of-fluid method. The code requires the installation of two libraries - AMReX and IRL, 
both of which are provided in this repository. Hence, this repository is 
self-contained, and no external packages need to be installed.

# Water jet in supersonic air crossflow 
<img src="Images/LJSCF_AMReX.gif?raw=true&v=100" alt="your_alternative_text" width="50%" height="50%" loop="true" autoplay="true"><img src="Images/SprayAcoustic.gif?raw=true&v=100" alt="your_alternative_text" width="50%" height="50%" loop="true" autoplay="true">  

<img src="Images/LJSCF_AMReX_Blender.gif?raw=true&v=100" alt="your_alternative_text" width="50%" height="50%" loop="true" autoplay="true">


# Spray atomization with acoustic excitation
 

# Installation and compilation 
```
git clone https://github.com/nataraj2/MultiphaseSolverAMReX.git
cd MultiphaseSolverAMReX
sh install.sh
```
```
vi amrex/Tools/GNUMake/Make.local
```
In the first 5 lines of `Make.local`, change the path to `IRLDIR`, and the compilers `CXX`, `CC`, 
`FC`, `F90`. Make sure to use MPI compilers.

# Running
The following is to run the case of Zalesak's disk.
```
cd amrex/Tutorials/Amr/MultiphaseSolver_Zalesak_Elvira/Exec/Zalesak_Elvira
make -j
sh run_3d.sh
```
## Input file description
1. `geometry.prob_lo` and `geometry.prob_hi` define the domain
2. `amr.max_level` is the number of refinement levels above the base level. It is set to 2. 
Hence there will be a total of 3 levels (including base level).
3. `amr.regrid_int` is the frequency of regridding. Has to be 1. Regridding has to be 
done at every step since we do not allow for the interface to cut across levels. The interface 
should always be enclosed by the finest level.
4. `amr.plot_int` and `amr.chk_int` are the frequencies at which the output 
and checkpoint (restart) files are written.
5. `amr.blocking_factor` is the blocking factor. Defines the size of the blocks.

# Visualization
VisIt 2.13 is the minimum version needed for visualization. The earlier versions of VisIt
cannot read the AMReX output files. In the running directory, do
```
sh run_output.sh
```
This script will create a text file `movie.visit` (which contains a list of headers of the solution files - which has all 
the flow variables). This also creates a directory `ensight-3D`, in which all plic data is 
written. A file `plic.visit` is also written inside `ensight-3D`. Both `movie.visit` and `plic.visit`
 can be loaded simultaneously into VisIt and visualized. The steps are as follows

1. `sh run_output.sh`
2. File-> Open file -> movie.visit
3. File-> Open file -> enisght-3D/plic.visit
4. Active source -> movie.visit. Add->Mesh->mesh
5. Active source -> plic.visit. Add->Mesh->mesh
6. VisIt will ask if a correlation should be setup. Click Yes. This will allow the plic and the mesh to 
change simultaneously
7. Active time slider -> Correlation01
8. Click Draw
9. Click the right arrow button below the Active time slider to animate




