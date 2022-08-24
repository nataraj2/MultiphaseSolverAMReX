module helmholtz

  use amrex_base_module
  use amrex_linear_solver_module

  implicit none

  ! parameters

  !integer, save :: max_level = 1
  !integer, save :: ref_ratio = 2
  !integer, save :: n_cell = 128
  !integer, save :: max_grid_size = 64

  logical, save :: composite_solve = .true.
  !logical, save :: composite_solve = .false.

  ! prob_type 1 here is Poisson with homogeneous Dirichlet boundary.
  ! prob_type 2 here is ABecLaplacian with homogeneous Neumann boundary.
  !integer, save :: prob_type = 1

  integer, save :: verbose = 2
  integer, save :: cg_verbose = 0
  integer, save :: max_iter = 100
  integer, save :: max_fmg_iter = 0
  integer, save :: linop_maxorder = 2
  logical, save :: agglomeration = .true.
  logical, save :: consolidation = .true.

  ! data
  !type(amrex_boxarray), allocatable, save :: basolve(:)
  !type(amrex_distromap), allocatable, save :: dmsolve(:)

  type(amrex_multifab), allocatable, save :: solution(:)
  type(amrex_multifab), allocatable, save :: rhs(:)
  type(amrex_multifab), allocatable, save :: exact_solution(:)
  type(amrex_multifab), allocatable, save :: acoef(:)
  type(amrex_multifab), allocatable, save :: bcoef(:)
  type(amrex_multifab),allocatable :: beta(:,:)
  real(amrex_real), save :: ascalar, bscalar
  type(amrex_abeclaplacian) :: abeclap
  type(amrex_multifab), allocatable, save :: gradP(:,:)

end module helmholtz
