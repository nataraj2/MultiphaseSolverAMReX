
module amr_data_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_base_module
  use bc_module
  use amrex_linear_solver_module 
  use amrex_fi_mpi
  use irl_fortran_interface

  implicit none

  logical, save :: composite_solve = .true.

  integer, save :: verbose = 2
  integer, save :: cg_verbose = 0
  integer, save :: max_iter = 10000
  integer, save :: max_fmg_iter = 0
  integer, save :: linop_maxorder = 2
  logical, save :: agglomeration = .true.
  logical, save :: consolidation = .true.

  double precision, parameter :: VOFlo=0.0000000001d0
  double precision, parameter :: VOFhi=0.9999999999d0

  real(rt), allocatable :: t_new(:)
  real(rt), allocatable :: t_old(:)
  real(amrex_real) :: dt_step

  type(amrex_multifab), allocatable :: phi_new(:)
  type(amrex_multifab), allocatable :: phi_old(:)
  type(amrex_multifab), allocatable :: phi_n(:)
  type(amrex_multifab), allocatable :: phi_mp(:)
  type(amrex_multifab), allocatable :: phi_mp_n(:)
  type(amrex_multifab), allocatable :: gradxphi(:)
  type(amrex_multifab), allocatable :: gradyphi(:)
  type(amrex_multifab), allocatable :: gradzphi(:)
  type(amrex_multifab), allocatable :: checkvar(:)
  type(amrex_multifab), allocatable :: norm(:)
  type(amrex_multifab), allocatable :: Gvolmfab(:),Lvolmfab(:)

  type tile_type
        real(amrex_real), allocatable, dimension(:) :: x,y,z,xm,ym,zm
  end type tile_type
  type(tile_type), target, allocatable :: tile(:)  
  type(amrex_multifab), allocatable :: grid_new(:)

  type(amrex_fluxregister), allocatable :: flux_reg(:)

  type(amrex_multifab), allocatable, save :: Pmfab(:)
  type(amrex_multifab), allocatable, save :: rhs(:)
  type(amrex_multifab), allocatable, save :: exact_solution(:)
  type(amrex_multifab), allocatable, save :: acoef(:)
  type(amrex_multifab), allocatable, save :: bcoef(:)
  type(amrex_multifab),allocatable :: beta(:,:), betatemp(:,:)
  real(amrex_real), save :: ascalar, bscalar
  type(amrex_abeclaplacian) :: abeclap
  type(amrex_multifab), allocatable, save :: gradPmfab(:,:), facevelmfab(:,:), Pfacemfab(:,:), rhofacemfab(:,:)
  type(amrex_multifab), allocatable, save :: Umfab(:),Vmfab(:),Wmfab(:),vortmagmfab(:) 
  type(CapDod_type), dimension(1) :: flux_hexahedron
  type(TagAccVM_SepVM_type) :: f_moments

  integer :: ncomp

  integer, allocatable, save :: stepno_vec(:)
  real(rt), allocatable, save :: dt_vec(:)
  logical :: do_reflux  = .true.

  type test_type
	double precision :: xx,yy
  end type test_type

  integer :: myrank, ierr1

contains

  subroutine amr_data_init ()

    allocate(t_new(0:amrex_max_level))
    t_new = 0.0_rt

    allocate(t_old(0:amrex_max_level))
    t_old = -1.0e100_rt

    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))
    allocate(phi_n(0:amrex_max_level))
    allocate(phi_mp(0:amrex_max_level))
    allocate(phi_mp_n(0:amrex_max_level))
    allocate(gradxphi(0:amrex_max_level))
    allocate(gradyphi(0:amrex_max_level))
    allocate(gradzphi(0:amrex_max_level))
    allocate(checkvar(0:amrex_max_level))
    allocate(flux_reg(0:amrex_max_level))
    allocate(grid_new(0:amrex_max_level))
    allocate(norm(0:amrex_max_level)) 
    allocate(Gvolmfab(0:amrex_max_level)) 
    allocate(Lvolmfab(0:amrex_max_level)) 

    allocate(Pmfab(0:amrex_max_level))
    allocate(rhs(0:amrex_max_level))
    allocate(exact_solution(0:amrex_max_level))
    !if (prob_type .eq. 2) then
    allocate(acoef(0:amrex_max_level))
    allocate(bcoef(0:amrex_max_level))

    allocate(beta(amrex_spacedim,0:amrex_max_level))
    allocate(betatemp(amrex_spacedim,0:amrex_max_level))
    allocate(gradPmfab(amrex_spacedim,0:amrex_max_level))
    allocate(facevelmfab(amrex_spacedim,0:amrex_max_level))
    allocate(Pfacemfab(amrex_spacedim,0:amrex_max_level))
    allocate(rhofacemfab(amrex_spacedim,0:amrex_max_level))

    allocate(Umfab(0:amrex_max_level))
    allocate(Vmfab(0:amrex_max_level))
    allocate(Wmfab(0:amrex_max_level))
    allocate(vortmagmfab(0:amrex_max_level))

    allocate(stepno_vec(0:amrex_max_level))
    allocate(dt_vec(0:amrex_max_level))

   call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr1)

  end subroutine amr_data_init

  subroutine amr_data_finalize
    integer :: lev
    do lev = 0, amrex_max_level
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       call amrex_multifab_destroy(phi_n(lev))
       call amrex_multifab_destroy(gradxphi(lev))
       call amrex_multifab_destroy(gradyphi(lev))
       call amrex_multifab_destroy(gradzphi(lev))
    end do
    do lev = 1, amrex_max_level
       call amrex_fluxregister_destroy(flux_reg(lev))
    end do
  end subroutine amr_data_finalize
  
end module amr_data_module
