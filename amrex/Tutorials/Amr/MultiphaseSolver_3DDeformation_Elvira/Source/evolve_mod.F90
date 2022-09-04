module evolve_module

  use amrex_amr_module
  use solve_helmholtz
  use rhs_helmholtz
  use mod_compute_cell_value
  use mod_compute_correction_step
  use mod_compute_interp_pressure
  use mod_compute_interp_rho
  use mod_compute_multiphase_variables
  use mod_compute_plic_half
  use mod_compute_vortmag

  implicit none
  private

  public :: evolve

contains

  subroutine evolve ()
    use my_amr_module, only : stepno, max_step, stop_time, dt, plot_int, chk_int, restart
    use amr_data_module
    use compute_dt_module, only : compute_dt
    use face_velocity_multifab
    use plotfile_module, only : writeplotfile, write_plotfile_multi2
    use restart_module
    use plic_module
    use irl_fortran_interface
  
    real(amrex_real) :: cur_time
    integer :: last_plot_file_step, step, lev, substep, finest_level
    integer :: icomp
    integer :: idim, nlevs
    integer :: cellcount(1), totalcellcount

    if(len_trim(restart).ne.0)then
      do lev=0,amrex_get_finest_level()
         stepno(lev)=stepno_vec(lev)
         dt(lev)=dt_vec(lev)
      enddo
    endif

    cur_time = t_new(0)
    last_plot_file_step = 0;
   
    call write_plotfile_multi2()
    !call compute_helmholtz_rhs
    call compute_VOF2Mesh(0)

    do step = stepno(0), max_step-1
       if (cur_time .ge. stop_time) exit

       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       call compute_dt()

       lev = 0
       substep = 1
       call timestep(lev, cur_time, substep)

       !totalcellcount=0
       !do lev=0, amrex_get_finest_level()
	!  cellcount=0	
        !  call countcellsatlevel(phi_new(lev)%ba%p,cellcount)
	!  totalcellcount=totalcellcount+cellcount(1)
       !end do
       ! if (amrex_parallel_ioprocessor()) then
       !open(unit=10,file='total_no_cells.dat',position='append')
	!write(10,*) step, totalcellcount
	!close(10)
	!endif
 
       cur_time = cur_time + dt(0)

       if (amrex_parallel_ioprocessor()) then
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", dt(0)
       end if

	dt_step=dt(0)

        finest_level = amrex_get_finest_level()
	
       ! sync up time to avoid roundoff errors
       finest_level = amrex_get_finest_level()
       do lev = 0, finest_level
          t_new(lev) = cur_time
       end do

	if(stepno(0).eq.1)then
	!call write_plotfile_multi2()
        !call compute_VOF2Mesh(stepno(0))
	endif


       if (plot_int .gt. 0 .and. mod(step+1,plot_int) .eq. 0) then
          last_plot_file_step = step+1
          call write_plotfile_multi2()
          call compute_VOF2Mesh(stepno(0))
       end if

       if(chk_int .gt. 0 .and. mod(step+1,chk_int) .eq. 0)then
       call writecheckpointfile (stepno,finest_level,dt,t_new,phi_new(0:finest_level)%ba%p,phi_new(0:finest_level)%p)
       endif

       if (cur_time .ge. stop_time - 1.e-6_amrex_real*dt(0)) exit

       nlevs=amrex_get_finest_level()

    end do

    if (plot_int .gt. 0 .and. stepno(0) .gt. last_plot_file_step) then
       call write_plotfile_multi2()
       call compute_VOF2Mesh(stepno(0))
    end if

  end subroutine evolve

  recursive subroutine timestep (lev, time, substep)
    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, flux_reg, do_reflux
    use averagedown_module, only : averagedownto
    use restart_module
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time
    integer :: cellcount(1)
    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep

    if (regrid_int .gt. 0) then
       if (.not.allocated(last_regrid_step)) then
          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
       end if

       ! regrid doesn't change the base level, so we don't regrid on amrex_max_level
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then
          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time) ! the finest level may change during regrid
             finest_level = amrex_get_finest_level()

             do k = lev, finest_level
                last_regrid_step(k) = stepno(k)
             end do

             do k = old_finest_level+1, finest_level
                dt(k) = dt(k-1) / amrex_ref_ratio(k-1)             
             end do
          end if
       end if
    end if

    stepno(lev) = stepno(lev)+1

    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpath.
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    ! swap phi_new(lev) and phi_old(lev) so they are consistent with t_new(lev) and t_old(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))

    call advance(lev, time, dt(lev), stepno(lev), substep, nsubsteps(lev))
   
    !call countcellsatlevel(phi_new(lev)%ba%p,cellcount)

    if (lev .lt. amrex_get_finest_level()) then
       do fine_substep = 1, nsubsteps(lev+1)
          call timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep)
       end do

       if (do_reflux) then
          call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real)
       end if

       call averagedownto(lev)
    end if    
  end subroutine timestep

  ! update phi_new(lev)
  subroutine advance (lev, time, dt, step, substep, nsub)
    use my_amr_module, only : verbose
    use amr_data_module
    use face_velocity_module
    use advect_module, only : advect
    use fillpatch_module, only : fillpatch
    use reconstruct_module
    use grid_module
    use irl_fortran_interface
    integer, intent(in) :: lev, step, substep, nsub
    real(amrex_real), intent(in) :: time, dt

    integer, parameter :: ngrow = 3
    integer :: idim
    logical :: nodal(3)
    type(amrex_multifab) :: phiborder,Uborder,phi_nborder,gradxphiborder,gradyphiborder,gradzphiborder,Gvolmfab_lev,Lvolmfab_lev
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx, tbx, bxnodal
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout,pux,puy,puz,pfx,pfy,pfz, &
         pf, pfab, phi_nborptr, pgradxphi,pgradyphi,pgradzphi,pnorm,Gvolptr,Lvolptr
    type(amrex_fab) :: facevel(amrex_spacedim)
    type(amrex_fab) ::  flux(amrex_spacedim)
    type(amrex_multifab) :: fluxes(amrex_spacedim)

    integer :: totaltilecount, tilecount
    integer :: ind(4),flo(3),fhi(3),clo(3),chi(3),normlo(4),normhi(4)
    type(ObjServer_PlanarLoc_type) :: planar_localizer_allocation
    type(ObjServer_PlanarSep_type) :: planar_separator_allocation
    type(ObjServer_LocSepLink_type) :: localized_separator_link_allocation
    type(ObjServer_LocLink_type) :: localizer_link_allocation
    type(PlanarLoc_type), allocatable, dimension(:,:,:) :: plocalizer
    !type(planarlocalizer_type), contiguous, pointer, dimension(:,:,:) :: plocalizer
    type(PlanarSep_type), allocatable, dimension(:,:,:), target :: pseparator
    type(LocSepLink_type), allocatable, dimension(:,:,:) :: plocseplink
    type(LocLink_type), dimension(:,:,:), allocatable :: ploclink
    real(amrex_real), contiguous, pointer, dimension(:) :: x,y,z,xm,ym,zm
    integer(IRL_LargeOffsetIndex_t) :: heapsize
    integer :: i,j,k
    integer :: allocateint

    if (verbose .gt. 0 .and. amrex_parallel_ioprocessor()) then
       write(*,'(A, 1X, I0, 1X, A, 1X, I0, A, 1X, G0)') &
            "[Level", lev, "step", step, "] ADVANCE with dt =", dt
    end if

    ncomp = phi_new(lev)%ncomp()

    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do
    end if

    !call amrex_multifab_build(grid_new(lev), phi_new(lev)%ba, phi_new(lev)%dm, 0, ngrow) 
    
  ! Compute the grid for all the tiles in this level
    call compute_grid(lev)

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    call amrex_multifab_build(phi_nborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)

    call fillpatch(lev, time, phiborder,phi_new,phi_old)
    call fillpatch(lev, time, phi_nborder,phi_n,phi_n) 
 
    ! Compute the reconstruction and gradients for update. The gradients are carried through the adaptation like phi_new
    ! and hence defined in amr_data_module
    call compute_reconstruct(lev,phiborder)

    ! Create and fill the multifabs for the gradients
    call amrex_multifab_build(gradxphiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    call amrex_multifab_build(gradyphiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    call amrex_multifab_build(gradzphiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)
    call fillpatch(lev, time, gradxphiborder, gradxphi, gradxphi)
    call fillpatch(lev, time, gradyphiborder, gradyphi, gradyphi)
    call fillpatch(lev, time, gradzphiborder, gradzphi, gradzphi)

    call amrex_multifab_build(Gvolmfab_lev, phi_new(lev)%ba, phi_new(lev)%dm, 6, 1)
    call amrex_multifab_build(Lvolmfab_lev, phi_new(lev)%ba, phi_new(lev)%dm, 6, 1)

    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    tilecount=0
    do while(mfi%next())
       tilecount=tilecount+1
       bx = mfi%tilebox()
       bxnodal = mfi%nodaltilebox()

       pin  => phiborder%dataptr(mfi)
       pout => phi_new(lev)%dataptr(mfi)
        
       phi_nborptr  => phi_nborder%dataptr(mfi)

       pgradxphi => gradxphiborder%dataptr(mfi)
       pgradyphi => gradyphiborder%dataptr(mfi)
       pgradzphi => gradzphiborder%dataptr(mfi)

       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call flux(idim)%resize(tbx,ncomp)
          call tbx%grow(1)
          call facevel(idim)%resize(tbx,1)
       end do

       pux => facevel(1)%dataptr()
       pfx =>  flux(1)%dataptr()
       puy => facevel(2)%dataptr()
       pfy =>  flux(2)%dataptr()
       puz => facevel(3)%dataptr()
       pfz =>  flux(3)%dataptr()

       ind=lbound(pin);clo(1)=ind(1);clo(2)=ind(2);clo(3)=ind(3)
       ind=ubound(pin);chi(1)=ind(1);chi(2)=ind(2);chi(3)=ind(3)

        x  => tile(tilecount)%x
        y  => tile(tilecount)%y
        z  => tile(tilecount)%z
        xm => tile(tilecount)%xm
        ym => tile(tilecount)%ym
        zm => tile(tilecount)%zm


       ! Initialize size for IRL
       heapsize = int((bx%hi(3)-bx%lo(3)+3),8)*int((bx%hi(2)-bx%lo(2)+3),8)*int((bx%hi(1)-bx%lo(1)+3),8)
       call new(planar_localizer_allocation,heapsize)
       call new(planar_separator_allocation,heapsize)
       call new(localized_separator_link_allocation,heapsize)
       call new(localizer_link_allocation,heapsize)

  	allocate(plocalizer(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)))
  	allocate(pseparator(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)))
  	allocate(plocseplink(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)))
  	allocate(ploclink(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)))

        !plocalizer = grid_new(lev)%localizer(mfi)
        !pseparator = grid_new(lev)%separator(mfi)
        !plocseplink = grid_new(lev)%locseplink(mfi)

	

        !if(step.eq.1)then
        call initialize_IRL_arrays(bx%lo,bx%hi,planar_localizer_allocation,planar_separator_allocation,localized_separator_link_allocation,localizer_link_allocation,&
				   plocalizer,pseparator,plocseplink,ploclink,clo,chi,x,y,z)

	Gvolptr => Gvolmfab_lev%dataptr(mfi)
        Lvolptr => Lvolmfab_lev%dataptr(mfi)
	pnorm => norm(lev)%dataptr(mfi)
        normlo=lbound(pnorm)
        normhi=ubound(pnorm)
        !if(stepno(lev).gt.1)then
        !call vof_normal(bx%lo,bx%hi,pin(:,:,:,ncomp),clo,chi,pseparator,normlo,normhi,pnorm,Gvolptr,Lvolptr,lbound(Gvolptr),ubound(Gvolptr),x,y,z,amrex_geom(lev)%dx)
        call vof_normal(bx%lo,bx%hi,pin(:,:,:,ncomp),clo,chi,pseparator,normlo,normhi,pnorm,Gvolptr,Lvolptr,lbound(Gvolptr),ubound(Gvolptr),amrex_geom(lev)%dx,amrex_problo)

	call get_face_velocity(phi_nborptr(:,:,:,1),phi_nborptr(:,:,:,2),phi_nborptr(:,:,:,3),phi_nborptr(:,:,:,4),&
                               phi_nborptr(:,:,:,6),phi_nborptr(:,:,:,7),phi_nborptr(:,:,:,8),phi_nborptr(:,:,:,9),phi_nborptr(:,:,:,11),&
                              lbound(phi_nborptr),ubound(phi_nborptr), &
                              pux, lbound(pux), ubound(pux), &
                              puy, lbound(puy), ubound(puy), &
                              puz, lbound(puz), ubound(puz), &
                              Gvolptr, Lvolptr, lbound(Gvolptr), ubound(Gvolptr),&
                              amrex_geom(lev)%dx, amrex_problo,time)


	!print*, "here 1"
	!call checkplane(bx%lo,bx%hi,pseparator2,clo,chi)

       call advect(time, bx%lo, bx%hi, &
            pin, lbound(pin), ubound(pin), &
            pout,lbound(pout),ubound(pout), &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
            puz, lbound(puz), ubound(puz), &
            pfx, lbound(pfx), ubound(pfx), &
            pfy, lbound(pfy), ubound(pfy), &
            pfz, lbound(pfz), ubound(pfz), &
            pgradxphi, pgradyphi, pgradzphi, &
	    plocalizer,pseparator,plocseplink,ploclink,&
            x,y,z,xm,ym,zm,amrex_geom(lev)%dx, amrex_problo, dt)

	call ObjServer_PlanarLoc_class_delete(planar_localizer_allocation)
	call ObjServer_PlanarSep_class_delete(planar_separator_allocation)
	call ObjServer_LocSepLink_class_delete(localized_separator_link_allocation)
	call ObjServer_LocLink_class_delete(localizer_link_allocation)
		

	!print*, "here 2"
	deallocate(plocalizer)
	deallocate(pseparator)
	deallocate(plocseplink)
	deallocate(ploclink)
	
       if (do_reflux) then
          do idim = 1, amrex_spacedim
             pf => fluxes(idim)%dataptr(mfi)
             pfab => flux(idim)%dataptr()
             tbx = mfi%nodaltilebox(idim)
             pf       (tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:) = &
                  pfab(tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:)
          end do
       end if
    end do
 
    call amrex_mfiter_destroy(mfi)
  
    deallocate(tile)

    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(facevel(idim))
       call amrex_fab_destroy( flux(idim))
    end do
    !$omp end parallel

    if (do_reflux) then
       ! Note that the fluxes have already been scaled by dt and area.
       if (lev > 0) then
          call flux_reg(lev)%fineadd(fluxes, 1.0_amrex_real)
       end if

       if (lev < amrex_get_finest_level()) then
          call flux_reg(lev+1)%crseinit(fluxes, -1.0_amrex_real)
       end if

       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(fluxes(idim))
       end do
    end if


    call amrex_multifab_destroy(phiborder)
    call amrex_multifab_destroy(phi_nborder)
    call amrex_multifab_destroy(gradxphiborder)
    call amrex_multifab_destroy(gradyphiborder)
    call amrex_multifab_destroy(gradzphiborder)
    call amrex_multifab_destroy(Gvolmfab_lev)
    call amrex_multifab_destroy(Lvolmfab_lev)
    !call amrex_multifab_destroy(grid_new(lev))

  end subroutine advance

end module evolve_module
