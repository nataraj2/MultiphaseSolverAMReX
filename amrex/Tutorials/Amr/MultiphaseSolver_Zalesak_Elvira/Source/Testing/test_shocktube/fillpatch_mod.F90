module fillpatch_module

  use iso_c_binding
  use amrex_amr_module
  use amr_data_module
  use bc_module
  use my_amrex_filcc_module

  implicit none

  private

  public :: fillpatch, fillcoarsepatch, fillpatchsingle

contains

  ! Fill phi with data from Grho_old and Grho of current level and one level below.
  subroutine fillpatch (lev, time, varmfab_border, varmfab, varmfab_old)
    use amr_data_module
    use bc_module, only : lo_bc, hi_bc
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: varmfab_border, varmfab(0:), varmfab_old(0:)
    
    !integer, parameter :: src_comp=2, dst_comp=2, num_comp=1  ! for this test code
    integer, parameter :: num_comp=1  ! for this test code
    integer :: src_comp, dst_comp

    do src_comp=1,ncomp
       dst_comp=src_comp       

    if (lev .eq. 0) then
       call amrex_fillpatch(varmfab_border, t_old(lev), varmfab_old(lev), &
            &                    t_new(lev), varmfab(lev), &
            &               amrex_geom(lev), fill_physbc , &
            &               time, src_comp, dst_comp, num_comp)
    else
       call amrex_fillpatch(varmfab_border, t_old(lev-1), varmfab_old(lev-1), &
            &                    t_new(lev-1), varmfab(lev-1), &
            &               amrex_geom(lev-1), fill_physbc   , &
            &                    t_old(lev  ), varmfab_old(lev  ), &
            &                    t_new(lev  ), varmfab(lev  ), &
            &               amrex_geom(lev  ), fill_physbc   , &
            &               time, src_comp, dst_comp, num_comp, &
            &               amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
            &               lo_bc, hi_bc)
       ! see amrex_interpolater_module for a list of interpolaters
    end if

    enddo
  end subroutine fillpatch

  subroutine fillcoarsepatch (lev, time, varmfab_border, varmfab, varmfab_old)
    use amr_data_module
    use bc_module, only : lo_bc, hi_bc
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: varmfab_border, varmfab(0:), varmfab_old(0:)

    !integer, parameter :: src_comp=2, dst_comp=2, num_comp=1  ! for this test code
    integer, parameter :: num_comp=1  ! for this test code
    integer :: src_comp, dst_comp

     print*, "Inside coarse patch"

    do src_comp=1,ncomp
       dst_comp=src_comp       

    call amrex_fillcoarsepatch(varmfab_border, t_old(lev-1), varmfab_old(lev-1),  &
         &                          t_new(lev-1), varmfab(lev-1),  &
         &                     amrex_geom(lev-1),    fill_physbc,  &
         &                     amrex_geom(lev  ),    fill_physbc,  &
         &                     time, src_comp, dst_comp, num_comp, &
         &                     amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
         &                     lo_bc, hi_bc)
       ! see amrex_interpolater_module for a list of interpolaters

     enddo

  end subroutine fillcoarsepatch

  subroutine fill_physbc (pmf, scomp, ncomp, time, pgeom) bind(c)
    use amrex_geometry_module, only : amrex_is_all_periodic
    use amrex_filcc_module, only : amrex_filcc
    use bc_module, only : lo_bc, hi_bc
    use my_amrex_filcc_module, only : my_amrex_filcc
    type(c_ptr), value :: pmf, pgeom
    integer(c_int), value :: scomp, ncomp
    real(amrex_real), value :: time

    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
    integer :: plo(4), phi(4)

    if (.not. amrex_is_all_periodic()) then

       geom = pgeom
       mf = pmf

       !$omp parallel private(mfi,p,plo,phi)
       call amrex_mfiter_build(mfi, mf, tiling=.false.)
       do while(mfi%next())
          p => mf%dataptr(mfi)
          if (.not. geom%domain%contains(p)) then ! part of this box is outside the domain
             plo = lbound(p)
             phi = ubound(p)
             call amrex_filcc(p, plo, phi,         & ! fortran array and bounds
                  geom%domain%lo, geom%domain%hi,  & ! index extent of whole problem domain
                  geom%dx,                         & ! cell size in real
                  geom%get_physical_location(plo), & ! physical location of lower left corner
                  lo_bc, hi_bc)                      ! bc types for each component

	     !call my_amrex_filcc(p, plo, phi,         & ! fortran array and bounds
              !    geom%domain%lo, geom%domain%hi,  & ! index extent of whole problem domain
               !   geom%dx,                         & ! cell size in real
                !  geom%get_physical_location(plo), & ! physical location of lower left corner
                 ! lo_bc, hi_bc)                      ! bc types for each component

             ! amrex_filcc doesn't fill EXT_DIR (see amrex_bc_types_module for a list of bc types
             ! In that case, the user needs to fill it.
          end if
       end do
       !$omp end parallel
       call amrex_mfiter_destroy(mfi)

    end if
 
  end subroutine fill_physbc

subroutine fillpatchsingle (lev, time, varmfab_border, varmfab, varmfab_old)
    use amr_data_module
    use bc_module, only : lo_bc, hi_bc
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: varmfab_border, varmfab(0:), varmfab_old(0:)
    
    !integer, parameter :: src_comp=2, dst_comp=2, num_comp=1  ! for this test code
    integer, parameter :: num_comp=1  ! for this test code
    integer :: src_comp, dst_comp

    do src_comp=1,1
       dst_comp=src_comp       

    if (lev .eq. 0) then
       call amrex_fillpatch(varmfab_border, t_old(lev), varmfab_old(lev), &
            &                    t_new(lev), varmfab(lev), &
            &               amrex_geom(lev), fill_physbc , &
            &               time, src_comp, dst_comp, num_comp)
    else
       call amrex_fillpatch(varmfab_border, t_old(lev-1), varmfab_old(lev-1), &
            &                    t_new(lev-1), varmfab(lev-1), &
            &               amrex_geom(lev-1), fill_physbc   , &
            &                    t_old(lev  ), varmfab_old(lev  ), &
            &                    t_new(lev  ), varmfab(lev  ), &
            &               amrex_geom(lev  ), fill_physbc   , &
            &               time, src_comp, dst_comp, num_comp, &
            &               amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
            &               lo_bc, hi_bc)
       ! see amrex_interpolater_module for a list of interpolaters
    end if

    enddo
  end subroutine fillpatchsingle

end module fillpatch_module
