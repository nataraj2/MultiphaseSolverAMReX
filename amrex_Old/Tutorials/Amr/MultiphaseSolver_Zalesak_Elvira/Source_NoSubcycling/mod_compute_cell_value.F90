module mod_compute_cell_value

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module

  implicit none

contains

subroutine compute_cell_value(facemfab)

    implicit none

    integer :: ilev
    integer :: clo(4), chi(4), flo(4), fhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phi, faceptr, checkptr
    type(amrex_multifab) :: facemfab(amrex_spacedim,0:amrex_max_level)
    integer :: nlevs

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs
        call amrex_multifab_build(checkvar(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, 1, 0)
    enddo 	

    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
          phi    => phi_new(ilev)%dataptr(mfi)
 	  checkptr => checkvar(ilev)%dataptr(mfi)
	  faceptr  => facemfab(2,ilev)%dataptr(mfi)
          clo = lbound(phi)
          chi = ubound(phi)
	  flo = lbound(faceptr)
	  fhi = ubound(faceptr)

          call cell_value(ilev, bx%lo, bx%hi, phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3),lbound(phi), ubound(phi), faceptr(:,:,:,1),&
  	       lbound(faceptr),ubound(faceptr),checkptr(:,:,:,1),amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

    end do

    !do ilev = 0, nlevs
    !	call amrex_multifab_destroy(checkvar(ilev))
    !end do

end subroutine compute_cell_value

subroutine cell_value(level, lo, hi, &
     Grho, GrhoU, GrhoE, clo, chi, faceval,&
     flo, fhi, check, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), clo(3), chi(3), flo(3), fhi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: Grho,GrhoU,GrhoE,check
  double precision, dimension(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)) :: faceval

  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
	   check(i,j,1) = (faceval(i,j,1)+faceval(i,j+1,1))/2.0d0
        end do
     end do
  end do

end subroutine cell_value
 
end module mod_compute_cell_value
