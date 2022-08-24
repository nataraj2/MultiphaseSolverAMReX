module mod_compute_multiphase_variables

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module

  implicit none

contains

subroutine build_multiphase_mfab(phi_mpmfab)

    use amr_data_module
    type(amrex_multifab) :: phi_mpmfab(0:amrex_max_level)
    integer :: lev, nlevs

    nlevs = amrex_get_finest_level()
	
    do lev=0, nlevs
       call amrex_multifab_build(phi_mpmfab(lev), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0)
    enddo	

end subroutine build_multiphase_mfab

subroutine compute_multiphase_variables(phimfab,phi_mpmfab)

    implicit none

    integer :: ilev
    integer :: flo(4), fhi(4), clo(4), chi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phi_mpptr, phiptr
    integer :: nlevs
    type(amrex_multifab) :: phimfab(0:amrex_max_level), phi_mpmfab(0:amrex_max_level)

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
	  phiptr   => phimfab(ilev)%dataptr(mfi)
	  phi_mpptr   => phi_mpmfab(ilev)%dataptr(mfi)

	  clo = lbound(phiptr)
	  chi = ubound(phiptr)

          call multiphase_variables(ilev, bx%lo, bx%hi, phiptr(:,:,:,1), phiptr(:,:,:,2), phiptr(:,:,:,3), phiptr(:,:,:,4), phiptr(:,:,:,5), &
							phiptr(:,:,:,6), phiptr(:,:,:,7), phiptr(:,:,:,8), phiptr(:,:,:,9), phiptr(:,:,:,10), phiptr(:,:,:,11), clo, chi, &
			                                phi_mpptr(:,:,:,1), phi_mpptr(:,:,:,2), phi_mpptr(:,:,:,3), phi_mpptr(:,:,:,4), phi_mpptr(:,:,:,5), &
							amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

    end do

end subroutine compute_multiphase_variables


subroutine multiphase_variables (level, lo, hi, Grho, GrhoU, GrhoV, GrhoW, GrhoE, Lrho, LrhoU, LrhoV, LrhoW, LrhoE, VOF, clo, chi, RHO, rhoU, rhoV, rhoW, rhoE, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), clo(3), chi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: Grho, GrhoU, GrhoV, GrhoW, GrhoE, Lrho, LrhoU, LrhoV, LrhoW, LrhoE, VOF
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: RHO, rhoU, rhoV, rhoW, rhoE
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
	   RHO(i,j,k)  = (1.0d0-VOF(i,j,k))*Grho(i,j,k)+VOF(i,j,k)*Lrho(i,j,k)
	   rhoU(i,j,k) = (1.0d0-VOF(i,j,k))*GrhoU(i,j,k)+VOF(i,j,k)*LrhoU(i,j,k)
	   rhoV(i,j,k) = (1.0d0-VOF(i,j,k))*GrhoV(i,j,k)+VOF(i,j,k)*LrhoV(i,j,k) 
	   rhoW(i,j,k) = (1.0d0-VOF(i,j,k))*GrhoW(i,j,k)+VOF(i,j,k)*LrhoW(i,j,k)
	   rhoE(i,j,k) = (1.0d0-VOF(i,j,k))*GrhoE(i,j,k)+VOF(i,j,k)*LrhoE(i,j,k)
        end do
     end do
  end do

end subroutine multiphase_variables

end module mod_compute_multiphase_variables
