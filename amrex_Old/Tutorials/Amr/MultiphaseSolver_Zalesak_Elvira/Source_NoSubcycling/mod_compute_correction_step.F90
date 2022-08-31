module mod_compute_correction_step

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module

  implicit none

contains

subroutine compute_correct_momentum

    implicit none

    integer :: ilev
    integer :: clo(4), chi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phiptr, PXptr, PYptr, PZptr
    integer :: nlevs

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
	  phiptr   => phi_mp(ilev)%dataptr(mfi)
          PXptr  => Pfacemfab(1,ilev)%dataptr(mfi)
          PYptr  => Pfacemfab(2,ilev)%dataptr(mfi)
          PZptr  => Pfacemfab(3,ilev)%dataptr(mfi)
	  clo = lbound(phiptr)
	  chi = ubound(phiptr)
          fxlo = lbound(PXptr)
          fxhi = ubound(PXptr)
          fylo = lbound(PYptr)
          fyhi = ubound(PYptr)
          fzlo = lbound(PZptr)
          fzhi = ubound(PZptr)

          call correct_momentum(ilev, bx%lo, bx%hi, phiptr(:,:,:,2), phiptr(:,:,:,3), phiptr(:,:,:,4), clo, chi, PXptr, fxlo, fxhi, &
				PYptr, fylo, fyhi, PZptr, fzlo, fzhi, amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

    end do

end subroutine compute_correct_momentum


subroutine correct_momentum (level, lo, hi, rhoU, rhoV, rhoW, clo, chi, PX, fxlo, fxhi, PY, fylo, fyhi, PZ, fzlo, fzhi, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3), clo(3), chi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: rhoU, rhoV, rhoW
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: PX
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: PY
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: PZ

  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
                rhoU(i,j,k) = rhoU(i,j,k) - dt_step/dx(1)*(PX(i+1,j,k)-PX(i,j,k))
                rhoV(i,j,k) = rhoV(i,j,k) - dt_step/dx(2)*(PY(i,j+1,k)-PY(i,j,k))
                rhoW(i,j,k) = rhoW(i,j,k) - dt_step/dx(3)*(PZ(i,j,k+1)-PZ(i,j,k))
		!GrhoU(i,j,k) = (PX(i,j,k)+PX(i+1,j,k))/2.0d0
        end do
     end do
  end do

end subroutine correct_momentum


subroutine compute_correct_face_velocities

    implicit none

    integer :: ilev
    integer :: fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: Ufaceptr, dPdxptr, RHOxptr, Vfaceptr, dPdyptr, RHOyptr, Wfaceptr, dPdzptr, RHOzptr
    integer :: nlevs

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
          Ufaceptr      => facevelmfab(1,ilev)%dataptr(mfi)
          Vfaceptr      => facevelmfab(2,ilev)%dataptr(mfi)
          Wfaceptr      => facevelmfab(3,ilev)%dataptr(mfi)
          dPdxptr  => gradPmfab(1,ilev)%dataptr(mfi)
          dPdyptr  => gradPmfab(2,ilev)%dataptr(mfi)
          dPdzptr  => gradPmfab(3,ilev)%dataptr(mfi)
	  RHOxptr   => rhofacemfab(1,ilev)%dataptr(mfi)
	  RHOyptr   => rhofacemfab(2,ilev)%dataptr(mfi)
	  RHOzptr   => rhofacemfab(3,ilev)%dataptr(mfi)
          fxlo = lbound(Ufaceptr)
          fxhi = ubound(Ufaceptr)
          fylo = lbound(Vfaceptr)
          fyhi = ubound(Vfaceptr)
          fzlo = lbound(Wfaceptr)
          fzhi = ubound(Wfaceptr)


          call correct_face_velocities(ilev, bx%lo, bx%hi, Ufaceptr, dPdxptr, RHOxptr, fxlo, fxhi, &
					                   Vfaceptr, dPdyptr, RHOyptr, fylo, fyhi, &
					                   Wfaceptr, dPdzptr, RHOzptr, fzlo, fzhi, &
							   amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

    end do

end subroutine compute_correct_face_velocities

subroutine correct_face_velocities (level, lo, hi, Uface, dPdx, RHOx, fxlo, fxhi, &
						   Vface, dPdy, RHOy, fylo, fyhi, &
						   Wface, dPdz, RHOz, fzlo, fzhi, &
						   dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: Uface, dPdx, RHOx
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: Vface, dPdy, RHOy
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: Wface, dPdz, RHOz

  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
                Uface(i,j,k) = Uface(i,j,k) - dt_step/dx(1)*dPdx(i,j,k)/RHOx(i,j,k)
        end do
     end do
  end do

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
                Vface(i,j,k) = Vface(i,j,k) - dt_step/dx(2)*dPdy(i,j,k)/RHOy(i,j,k)
        end do
     end do
  end do

  do k=lo(3),hi(3)+1
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
                Wface(i,j,k) = Wface(i,j,k) - dt_step/dx(3)*dPdz(i,j,k)/RHOz(i,j,k)
        end do
     end do
  end do

end subroutine correct_face_velocities

subroutine compute_correct_energy

    implicit none

    integer :: ilev
    integer :: clo(4), chi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phiptr, phi_mpptr, PXptr, Ufaceptr, PYptr, Vfaceptr, PZptr, Wfaceptr
    integer :: nlevs

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
	  phiptr    => phi_new(ilev)%dataptr(mfi)
	  phi_mpptr => phi_mp(ilev)%dataptr(mfi)
          PXptr     => Pfacemfab(1,ilev)%dataptr(mfi)
          PYptr     => Pfacemfab(2,ilev)%dataptr(mfi)
          PZptr     => Pfacemfab(3,ilev)%dataptr(mfi)
          Ufaceptr  => facevelmfab(1,ilev)%dataptr(mfi)
          Vfaceptr  => facevelmfab(2,ilev)%dataptr(mfi)
          Wfaceptr  => facevelmfab(3,ilev)%dataptr(mfi)
	  clo = lbound(phiptr)
	  chi = ubound(phiptr)
          fxlo = lbound(PXptr)
          fxhi = ubound(PXptr)
          fylo = lbound(PYptr)
          fyhi = ubound(PYptr)
          fzlo = lbound(PZptr)
          fzhi = ubound(PZptr)


          call correct_energy(ilev, bx%lo, bx%hi, phiptr(:,:,:,1), phiptr(:,:,:,5), phiptr(:,:,:,6), phiptr(:,:,:,10), phi_mpptr(:,:,:,1), clo, chi, PXptr, Ufaceptr, fxlo, fxhi, &
									     PYptr, Vfaceptr, fylo, fyhi, &
									     PZptr, Wfaceptr, fzlo, fzhi, &
									     amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

    end do

end subroutine compute_correct_energy

subroutine correct_energy (level, lo, hi, Grho, GrhoE, Lrho, LrhoE, RHO, clo, chi, PX, Uface, fxlo, fxhi,& 
							   PY, Vface, fylo, fyhi,&
							   PZ, Wface, fzlo, fzhi,&
							   dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), clo(3), chi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3)
  double precision, dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: Grho, GrhoE, Lrho, LrhoE, RHO
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: PX, Uface
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: PY, Vface
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: PZ, Wface

  double precision, intent(in) :: dx(3), prob_lo(3)
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
                GrhoE(i,j,k) = GrhoE(i,j,k)-Grho(i,j,k)*dt_step/RHO(i,j,k)*( (PX(i+1,j,k)*Uface(i+1,j,k)-PX(i,j,k)*Uface(i,j,k))/dx(1)+&
						      (PY(i,j+1,k)*Vface(i,j+1,k)-PY(i,j,k)*Vface(i,j,k))/dx(2)+&
						      (PZ(i,j,k+1)*Wface(i,j,k+1)-PZ(i,j,k)*Wface(i,j,k))/dx(3) )	
		LrhoE(i,j,k) = LrhoE(i,j,k)-Lrho(i,j,k)*dt_step/RHO(i,j,k)*( (PX(i+1,j,k)*Uface(i+1,j,k)-PX(i,j,k)*Uface(i,j,k))/dx(1)+&
						      (PY(i,j+1,k)*Vface(i,j+1,k)-PY(i,j,k)*Vface(i,j,k))/dx(2)+&
						      (PZ(i,j,k+1)*Wface(i,j,k+1)-PZ(i,j,k)*Wface(i,j,k))/dx(3) )	

                !GrhoE(i,j,k) = PX(i,j,k)
        end do
     end do
  end do

end subroutine correct_energy


subroutine compute_update_cons_momentum

    implicit none

    integer :: ilev
    integer :: clo(4), chi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: phiptr, phi_mpptr
    integer :: nlevs

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)

       do while (mfi%next())
          bx = mfi%tilebox()
          phi_mpptr => phi_mp(ilev)%dataptr(mfi)
          phiptr   => phi_new(ilev)%dataptr(mfi)
          clo = lbound(phiptr)
          chi = ubound(phiptr)

          call update_cons_momentum(ilev, bx%lo, bx%hi, phi_mpptr(:,:,:,1), phi_mpptr(:,:,:,2), phi_mpptr(:,:,:,3),phi_mpptr(:,:,:,4),clo, chi, &
				    phiptr(:,:,:,1), phiptr(:,:,:,2), phiptr(:,:,:,3), phiptr(:,:,:,4), &
				    phiptr(:,:,:,6), phiptr(:,:,:,7), phiptr(:,:,:,8), phiptr(:,:,:,9), &
			            amrex_geom(ilev)%dx, amrex_problo)

       end do

       call amrex_mfiter_destroy(mfi)

    end do

  end subroutine compute_update_cons_momentum

subroutine update_cons_momentum(level,lo,hi,RHO,rhoU,rhoV,rhoW,clo,chi,Grho,GrhoU,GrhoV,GrhoW,Lrho,LrhoU,LrhoV,LrhoW,dx,prob_lo)
   integer :: level, i, j, k
   integer, intent(in) :: lo(3), hi(3), clo(3), chi(3)
   double precision, intent(in), dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: RHO,rhoU,rhoV,rhoW
   double precision, intent(inout), dimension(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3)) :: Grho,GrhoU,GrhoV,GrhoW,Lrho,LrhoU,LrhoV,LrhoW
   real(amrex_real) :: dx(3), prob_lo(3)

   do k=clo(3), chi(3)
      do j=clo(2), chi(2)
         do i=clo(1), chi(1)
               GrhoU(i,j,k)=Grho(i,j,k)*rhoU(i,j,k)/RHO(i,j,k)
               GrhoV(i,j,k)=Grho(i,j,k)*rhoV(i,j,k)/RHO(i,j,k)
               GrhoW(i,j,k)=Grho(i,j,k)*rhoW(i,j,k)/RHO(i,j,k)
               LrhoU(i,j,k)=Lrho(i,j,k)*rhoU(i,j,k)/RHO(i,j,k)
               LrhoV(i,j,k)=Lrho(i,j,k)*rhoV(i,j,k)/RHO(i,j,k)
               LrhoW(i,j,k)=Lrho(i,j,k)*rhoW(i,j,k)/RHO(i,j,k)
         enddo
     enddo
  enddo

  end subroutine update_cons_momentum

end module mod_compute_correction_step
