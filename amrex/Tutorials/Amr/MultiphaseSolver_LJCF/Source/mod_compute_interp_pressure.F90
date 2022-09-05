module mod_compute_interp_pressure

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module
  use fillpatch_module

  implicit none

contains

subroutine compute_interp_pressure

    implicit none

    integer :: ilev
    integer :: prlo(4), prhi(4), borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4), vollo(4), volhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: PXptr, PYptr, PZptr, Pptr, phiborderptr, Gvolptr, Lvolptr
    integer :: nlevs
    integer, parameter :: ngrow = 3
    type(amrex_multifab) :: phiborder
    real(amrex_real) :: time
	
    time=0.0d0

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

       call amrex_multifab_build(phiborder, phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, ngrow)
       call fillpatch (ilev, time, phiborder, phi_new, phi_new)

       call amrex_mfiter_build(mfi, phi_new(ilev), tiling=.true.)
       do while (mfi%next())
          bx = mfi%tilebox()
 	  Pptr 	       => Pmfab(ilev)%dataptr(mfi)
 	  phiborderptr => phiborder%dataptr(mfi)
	  Gvolptr => Gvolmfab(ilev)%dataptr(mfi)
          Lvolptr => Lvolmfab(ilev)%dataptr(mfi)
	  PXptr      => Pfacemfab(1,ilev)%dataptr(mfi)
	  PYptr      => Pfacemfab(2,ilev)%dataptr(mfi)
	  PZptr      => Pfacemfab(3,ilev)%dataptr(mfi)
	  borlo  = lbound(phiborderptr)
	  borhi  = ubound(phiborderptr)
          prlo   = lbound(Pptr)
          prhi   = ubound(Pptr)
	  vollo  = lbound(Gvolptr)
	  volhi  = ubound(Gvolptr)
	  fxlo    = lbound(PXptr)
	  fxhi    = ubound(PXptr)
	  fylo    = lbound(PYptr)
	  fyhi    = ubound(PYptr)
	  fzlo    = lbound(PZptr)
	  fzhi    = ubound(PZptr)

        call interp_pressure(ilev, bx%lo, bx%hi, PXptr, fxlo,fxhi, &
			     PYptr,fylo,fyhi,PZptr,fzlo,fzhi,Pptr,prlo,prhi,&
	    	             phiborderptr(:,:,:,1),phiborderptr(:,:,:,6),borlo,borhi,Gvolptr,Lvolptr,vollo,volhi,&
	   	             amrex_geom(ilev)%dx, amrex_problo)

       end do
       call amrex_mfiter_destroy(mfi)
     
       call amrex_multifab_destroy(phiborder)

    end do

end subroutine compute_interp_pressure

subroutine interp_pressure(level, lo, hi, PX, fxlo, fxhi, &
			   PY, fylo, fyhi, PZ, fzlo, fzhi, P, prlo, prhi, &
			   Grho, Lrho, borlo, borhi, Gvol, Lvol, vollo, volhi, &
			   dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), prlo(3), prhi(3), borlo(3), borhi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3), vollo(4), volhi(4)
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: Grho, Lrho
  double precision, dimension(vollo(1):volhi(1),vollo(2):volhi(2),vollo(3):volhi(3),6) :: Gvol, Lvol
  double precision, dimension(prlo(1):prhi(1),prlo(2):prhi(2),prlo(3):prhi(3)) :: P
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: PX
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: PY
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: PZ
  double precision, intent(in) :: dx(3), prob_lo(3)
  double precision :: rho_l,rho_r,vol_l,vol_r
  integer          :: i,j,k
  double precision :: x,y,z

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
	   ! X face
           vol_r=Gvol(i,j,k,1)+Lvol(i,j,k,1); rho_r=Gvol(i,j,k,1)*Grho(i,j,k)+Lvol(i,j,k,1)*Lrho(i,j,k)
           vol_l=Gvol(i-1,j,k,2)+Lvol(i-1,j,k,2); rho_l=Gvol(i-1,j,k,2)*Grho(i-1,j,k)+Lvol(i-1,j,k,2)*Lrho(i-1,j,k)
	   PX(i,j,k)=(rho_r*P(i-1,j,k)+rho_l*P(i,j,k))/(rho_l+rho_r)
    	   !PX(i,j,k) = (P(i-1,j,k)*Grho(i,j,k)+P(i,j,k)*Grho(i-1,j,k))/(Grho(i-1,j,k)+Grho(i,j,k))
        end do
     end do
  end do

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
	   ! Y face
           vol_r=Gvol(i,j,k,3)+Lvol(i,j,k,3); rho_r=Gvol(i,j,k,3)*Grho(i,j,k)+Lvol(i,j,k,3)*Lrho(i,j  ,k)
           vol_l=Gvol(i,j-1,k,4)+Lvol(i,j-1,k,4); rho_l=Gvol(i,j-1,k,4)*Grho(i,j-1,k)+Lvol(i,j-1,k,4)*Lrho(i,j-1,k)
	   PY(i,j,k)=(rho_r*P(i,j-1,k)+rho_l*P(i,j,k))/(rho_l+rho_r)
    	   !PY(i,j,k) = (P(i,j-1,k)*Grho(i,j,k)+P(i,j,k)*Grho(i,j-1,k))/(Grho(i,j-1,k)+Grho(i,j,k))
        end do
     end do
  end do

  do k=lo(3),hi(3)+1
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
	   ! Z face
           vol_r=Gvol(i,j,k,5)+Lvol(i,j,k,5); rho_r=Gvol(i,j,k,5)*Grho(i,j,k)+Lvol(i,j,k,5)*Lrho(i,j,k)
           vol_l=Gvol(i,j,k-1,6)+Lvol(i,j,k-1,6); rho_l=Gvol(i,j,k-1,6)*Grho(i,j,k-1)+Lvol(i,j,k-1,6)*Lrho(i,j,k-1)
	   PZ(i,j,k)=(rho_r*P(i,j,k-1)+rho_l*P(i,j,k))/(rho_l+rho_r)
    	   !PZ(i,j,k) = (P(i,j,k-1)*Grho(i,j,k)+P(i,j,k)*Grho(i,j,k-1))/(Grho(i,j,k-1)+Grho(i,j,k))
        end do
     end do
  end do


end subroutine interp_pressure
 
end module mod_compute_interp_pressure
