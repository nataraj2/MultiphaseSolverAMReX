module my_amr_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real

  use amrex_base_module
  use amr_data_module

  use rhs_helmholtz

  implicit none

  ! runtime parameters
  !integer :: verbose    = 0

  integer :: max_step   = huge(1)
  integer :: regrid_int = 2
  integer :: chk_int  = -1
  integer :: plot_int   = -1
  !
  !logical :: do_reflux  = .true.
  !logical :: do_reflux  = .false.
  !
  real(rt) :: stop_time  = huge(1._rt)
  real(rt) :: cfl        = 0.7_rt
  !
  character(len=:), allocatable, save :: check_file
  character(len=:), allocatable, save :: plot_file
  character(len=:), allocatable, save :: restart

  integer, allocatable, save :: stepno(:)
  integer, allocatable, save :: nsubsteps(:)

  real(rt), allocatable, save :: dt(:)

  integer, private, parameter :: nghost = 0
 
contains

  subroutine my_amr_init ()
    use bc_module, only : lo_bc, hi_bc
    type(amrex_parmparse) :: pp
    integer :: ilev

    if (.not.amrex_amrcore_initialized()) call amrex_amrcore_init()
    
    call amrex_init_virtual_functions (my_make_new_level_from_scratch, &
         &                             my_make_new_level_from_coarse,  &
         &                             my_remake_level,                &
         &                             my_clear_level,                 &
         &                             my_error_estimate)

    ! some default parameters
    allocate(character(len=3)::check_file)
    check_file = "chk"
    allocate(character(len=3)::plot_file)
    plot_file = "plt"
    allocate(character(len=0)::restart)

    ! Read parameters
    call amrex_parmparse_build(pp)
    call pp%query("max_step", max_step)
    call pp%query("stop_time", stop_time)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters amr.*
    call amrex_parmparse_build(pp, "amr")
    call pp%query("regrid_int", regrid_int)
    call pp%query("ncomp", ncomp)
    call pp%query("chk_int", chk_int)
    call pp%query("plot_int", plot_int)
    call pp%query("check_file", check_file)
    call pp%query("plot_file", plot_file)
    call pp%query("restart", restart)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters myamr.*
    call amrex_parmparse_build(pp, "myamr")
    call pp%query("v", verbose)
    call pp%query("verbose", verbose)
    call pp%query("cfl", cfl)
    call pp%query("do_reflux", do_reflux)
    call amrex_parmparse_destroy(pp)

    
    allocate(stepno(0:amrex_max_level))
    stepno = 0

    allocate(nsubsteps(0:amrex_max_level))
    nsubsteps(0) = 1
    do ilev = 1, amrex_max_level
       nsubsteps(ilev) = amrex_ref_ratio(ilev-1)
    end do

    allocate(dt(0:amrex_max_level))
    dt = 1.d100

    allocate(lo_bc(amrex_spacedim,ncomp))
    allocate(hi_bc(amrex_spacedim,ncomp))

    !lo_bc(1:amrex_spacedim,1:ncomp) = amrex_bc_ext_dir
    !hi_bc(1:amrex_spacedim,1:ncomp) = amrex_bc_ext_dir

    !lo_bc = amrex_lo_periodic
    !hi_bc = amrex_lo_periodic
       lo_bc = amrex_bc_foextrap
       hi_bc = amrex_bc_foextrap 

       lo_bc(1,1:ncomp) = amrex_bc_ext_dir

    if (.not. amrex_is_all_periodic()) then
       !lo_bc = amrex_bc_foextrap
       !hi_bc = amrex_bc_foextrap 

       !lo_bc(1,:) = amrex_bc_foextrap ! left
        !lo_bc(1,:) = amrex_bc_ext_dir ! left
       !hi_bc(1,:) = amrex_bc_foextrap ! right
        !hi_bc(1,:) = amrex_bc_ext_dir ! right

       !lo_bc(2,:) = amrex_bc_foextrap ! top
        !lo_bc(2,:) = amrex_bc_ext_dir ! top
       !hi_bc(2,:) = amrex_bc_foextrap ! bottom
        !hi_bc(2,:) = amrex_bc_ext_dir ! bottom

        !lo_bc(3,:) = amrex_bc_ext_dir ! top
       !hi_bc(2,:) = amrex_bc_foextrap ! bottom
        !hi_bc(3,:) = amrex_bc_ext_dir ! bottom

    end if
	!lo_bc = amrex_bc_foextrap
       !hi_bc = amrex_bc_foextrap 

    call amr_data_init()

  end subroutine my_amr_init


  subroutine my_amr_finalize ()
    call amr_data_finalize()
  end subroutine my_amr_finalize

  ! Make a new level from scratch and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)
    use mod_compute_vortmag
    implicit none
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)
    integer :: idim
    logical :: nodal(3)
    integer :: icomp
    integer :: ngrow=3

    ba = pba
    dm = pdm

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call my_clear_level(lev)
  
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_n(lev), ba, dm, ncomp, nghost)

    call amrex_multifab_build(gradxphi(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(gradyphi(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(gradzphi(lev), ba, dm, ncomp, nghost)

    call amrex_multifab_build(grid_new(lev), ba, dm, 0, 3)
    call amrex_multifab_build(norm(lev), ba, dm, 3, 3)

   if (lev > 0 .and. do_reflux) then
      call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
   end if

    call init_prob_start(lev)
    
    do icomp=1,ncomp
       call phi_n(lev)%copy(phi_new(lev), icomp, icomp, 1, 0)
    enddo

  end subroutine my_make_new_level_from_scratch

  ! Make a new level from coarse level and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_coarse (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillcoarsepatch
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm

    integer :: idim
    logical :: nodal(3)

    type(amrex_multifab) :: null(0:amrex_max_level)
    integer :: icomp

    ba = pba
    dm = pdm

    call create_newlevel_and_copy(lev,time,phi_new,phi_old,ba,dm,ncomp,0,.false.,.true.)
    call create_newlevel_and_copy(lev,time,phi_n,phi_n,ba,dm,ncomp,0,.true.,.true.)
    call create_newlevel_and_copy(lev,time,gradxphi,gradxphi,ba,dm,ncomp,0,.true.,.true.)
    call create_newlevel_and_copy(lev,time,gradyphi,gradyphi,ba,dm,ncomp,0,.true.,.true.)
    call create_newlevel_and_copy(lev,time,gradzphi,gradzphi,ba,dm,ncomp,0,.true.,.true.)
    call create_newlevel_and_copy(lev,time,grid_new,grid_new,ba,dm,0,3,.true.,.false.)
    call create_newlevel_and_copy(lev,time,norm,norm,ba,dm,3,3,.true.,.false.)

    call my_clear_level(lev)

    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

  end subroutine my_make_new_level_from_coarse

  ! Remake a level from current and coarse elvels and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_remake_level (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillpatch
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: null(0:amrex_max_level)

    integer :: icomp
    integer :: idim
    logical :: nodal(3)

    ba = pba
    dm = pdm

    call create_remakelevel_and_copy(lev,time,phi_new,phi_old,ba,dm,ncomp,0,.false.,.true.)
    call create_remakelevel_and_copy(lev,time,phi_n,phi_n,ba,dm,ncomp,0,.true.,.true.)
    call create_remakelevel_and_copy(lev,time,gradxphi,gradxphi,ba,dm,ncomp,0,.true.,.true.)
    call create_remakelevel_and_copy(lev,time,gradyphi,gradyphi,ba,dm,ncomp,0,.true.,.true.)
    call create_remakelevel_and_copy(lev,time,gradzphi,gradzphi,ba,dm,ncomp,0,.true.,.true.)
    call create_remakelevel_and_copy(lev,time,grid_new,grid_new,ba,dm,0,3,.true.,.false.)
    call create_remakelevel_and_copy(lev,time,norm,norm,ba,dm,3,3,.true.,.false.)

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

  end subroutine my_remake_level

subroutine create_remakelevel_and_copy(lev,time,varmfab,varmfab_old,ba,dm,numcomp,ngrow,isidentical,isfillpatch)

    use fillpatch_module, only : fillpatch
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: new_phi_new, varmfab(0:), varmfab_old(0:)
    integer :: icomp
    integer :: idim
    logical, intent(in) :: isidentical,isfillpatch
    integer, intent(in) :: numcomp,ngrow

    call amrex_multifab_build(new_phi_new, ba, dm, numcomp, ngrow)
    if(isfillpatch.eqv..true.)then
    call fillpatch(lev, time, new_phi_new, varmfab, varmfab_old)
    endif

    call amrex_multifab_destroy(varmfab(lev))
    call amrex_multifab_build(varmfab(lev), ba, dm, numcomp, ngrow)
    if(isidentical.eqv..false.)then
    call amrex_multifab_destroy(varmfab_old(lev))
    call amrex_multifab_build(varmfab_old(lev), ba, dm, numcomp, ngrow)
    endif

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    do icomp=1,numcomp
       call varmfab(lev)%copy(new_phi_new, icomp, icomp, 1, ngrow)
    enddo

    call amrex_multifab_destroy(new_phi_new)

  end subroutine create_remakelevel_and_copy

subroutine create_newlevel_and_copy(lev,time,varmfab,varmfab_old,ba,dm,numcomp,ngrow,isidentical,isfillpatch)

    use fillpatch_module, only : fillcoarsepatch
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: new_phi_new, varmfab(0:), varmfab_old(0:)
    integer :: icomp
    integer :: idim
    logical, intent(in) :: isidentical,isfillpatch
    integer, intent(in) :: numcomp,ngrow

    call amrex_multifab_build(new_phi_new, ba, dm, numcomp, ngrow)
    if(isfillpatch.eqv..true.)then
    call fillcoarsepatch(lev, time, new_phi_new, varmfab, varmfab_old)
    endif
    call amrex_multifab_destroy(varmfab(lev))
    call amrex_multifab_build(varmfab(lev), ba, dm, numcomp, ngrow)
    if(isidentical.eqv..false.)then
    call amrex_multifab_destroy(varmfab_old(lev))
    call amrex_multifab_build(varmfab_old(lev), ba, dm, numcomp, ngrow)
    endif

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    do icomp=1,numcomp
       call varmfab(lev)%copy(new_phi_new, icomp, icomp, 1, ngrow)
    enddo

    call amrex_multifab_destroy(new_phi_new)

  end subroutine create_newlevel_and_copy

  subroutine my_clear_level (lev) bind(c)
    integer, intent(in), value :: lev
    integer :: idim
    call amrex_fluxregister_destroy(flux_reg(lev))

  end subroutine my_clear_level

  subroutine my_error_estimate (lev, cp, t, settag, cleartag) bind(c)
    use tagging_module, only : tag_phi_error
    use amr_data_module
    use mod_compute_vortmag
    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t
    character(kind=c_char), intent(in), value :: settag, cleartag

    real(amrex_real), allocatable, save :: phierr(:)
    type(amrex_parmparse) :: pp
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phiarr(:,:,:,:), vortmagptr(:,:,:,:)
    character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)

    if (.not.allocated(phierr)) then
       call amrex_parmparse_build(pp, "myamr")
       call pp%getarr("phierr", phierr)
       call amrex_parmparse_destroy(pp)
    end if

    tag = cp

    call amrex_multifab_build(vortmagmfab(lev), phi_new(lev)%ba, phi_new(lev)%dm, 1, 0)
    call compute_vortmag(lev)

    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       phiarr => phi_new(lev)%dataptr(mfi)
       vortmagptr => vortmagmfab(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
            phiarr, lbound(phiarr), ubound(phiarr), vortmagptr, &
            tagarr, lbound(tagarr), ubound(tagarr), &
            phierr(lev+1), settag, cleartag,amrex_geom(lev)%dx, amrex_problo)  ! +1 because level starts with 0, but phierr starts with 1
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

    call amrex_multifab_destroy(vortmagmfab(lev))

  end subroutine my_error_estimate
  
end module my_amr_module
