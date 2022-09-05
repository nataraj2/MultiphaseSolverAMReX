module mod_compute_liquid_mass

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module

  implicit none

contains

subroutine compute_liquid_mass

    implicit none

    integer :: ilev
    integer :: flo(4), fhi(4), clo(4), chi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phiptr
    integer :: nlevs
    double precision :: liq_mass, total_liq_mass

    nlevs = amrex_get_finest_level()

    total_liq_mass=0.0d0
    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
	  phiptr   => phi_new(ilev)%dataptr(mfi)

	  clo = lbound(phiptr)
	  chi = ubound(phiptr)

          call liquid_mass(ilev, bx%lo, bx%hi, phiptr(:,:,:,6),  phiptr(:,:,:,11), clo, chi, liq_mass, amrex_geom(ilev)%dx, amrex_problo)
		total_liq_mass=total_liq_mass+liq_mass

       end do

       call amrex_mfiter_destroy(mfi)

    end do
   
    print*,total_liq_mass

end subroutine compute_liquid_mass


subroutine liquid_mass (level, lo, hi, Lrho, VOF, clo, chi, liq_mass, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), clo(3), chi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: Lrho, VOF
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z, liq_mass

  liq_mass=0.0d0
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
		liq_mass=liq_mass+VOF(i,j,k)*dx(1)*dx(2)*dx(3)!*Lrho(i,j,k)
        end do
     end do
  end do

end subroutine liquid_mass

end module mod_compute_liquid_mass
