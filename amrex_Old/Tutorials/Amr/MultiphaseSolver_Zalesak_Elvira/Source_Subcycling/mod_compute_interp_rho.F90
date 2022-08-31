module mod_compute_interp_rho

  use amrex_base_module
  use amrex_amr_module
  use amr_data_module
  use fillpatch_module

  implicit none

contains

subroutine compute_interp_rho

    implicit none

    integer :: ilev
    integer :: borlo(4), borhi(4), fxlo(4), fxhi(4), fylo(4), fyhi(4), fzlo(4), fzhi(4), vollo(4), volhi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: RHOxptr, RHOyptr, RHOzptr, betaxptr, betayptr, betazptr, phiborderptr, Gvolptr, Lvolptr
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
 	  phiborderptr => phiborder%dataptr(mfi)
	  Gvolptr => Gvolmfab(ilev)%dataptr(mfi)
          Lvolptr => Lvolmfab(ilev)%dataptr(mfi)
	  RHOxptr  => rhofacemfab(1,ilev)%dataptr(mfi)
	  RHOyptr  => rhofacemfab(2,ilev)%dataptr(mfi)
	  RHOzptr  => rhofacemfab(3,ilev)%dataptr(mfi)
	  betaxptr  => beta(1,ilev)%dataptr(mfi)
	  betayptr  => beta(2,ilev)%dataptr(mfi)
	  betazptr  => beta(3,ilev)%dataptr(mfi)
          borlo = lbound(phiborderptr)
          borhi = ubound(phiborderptr)
	  vollo = lbound(Gvolptr)
          volhi = ubound(Gvolptr)
	  fxlo = lbound(RHOxptr)
	  fxhi = ubound(RHOxptr)
	  fylo = lbound(RHOyptr)
	  fyhi = ubound(RHOyptr)
	  fzlo = lbound(RHOzptr)
	  fzhi = ubound(RHOzptr)

	  call interp_rho(ilev, bx%lo, bx%hi, RHOxptr,betaxptr,&
               fxlo,fxhi,RHOyptr,betayptr,fylo,fyhi,RHOzptr,betazptr,fzlo,fzhi,phiborderptr(:,:,:,1),phiborderptr(:,:,:,6),borlo,borhi,Gvolptr,Lvolptr,vollo,volhi,amrex_geom(ilev)%dx, amrex_problo)

       end do
       call amrex_mfiter_destroy(mfi)
     
       call amrex_multifab_destroy(phiborder)

    end do

end subroutine compute_interp_rho

subroutine interp_rho(level, lo, hi, RHOx, betax, fxlo, fxhi, RHOy, betay, fylo, fyhi, RHOz, betaz, fzlo, fzhi, Grho, Lrho, borlo, borhi, Gvol, Lvol, vollo, volhi, dx, prob_lo)

  implicit none
  integer, intent(in) :: level, lo(3), hi(3), borlo(3), borhi(3), fxlo(3), fxhi(3), fylo(3), fyhi(3), fzlo(3), fzhi(3), vollo(4), volhi(4)
  double precision, dimension(borlo(1):borhi(1),borlo(2):borhi(2),borlo(3):borhi(3)) :: Grho, Lrho
  double precision, dimension(vollo(1):volhi(1),vollo(2):volhi(2),vollo(3):volhi(3),6) :: Gvol, Lvol
  double precision, dimension(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3)) :: RHOx, betax
  double precision, dimension(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3)) :: RHOy, betay
  double precision, dimension(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3)) :: RHOz, betaz
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
           RHOx (i,j,k)=(rho_l           +rho_r         )/(vol_l+vol_r)
	   betax(i,j,k)= 1.0d0/RHOx(i,j,k)
        end do
     end do
  end do

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
	   ! Y face
           vol_r=Gvol(i,j,k,3)+Lvol(i,j,k,3); rho_r=Gvol(i,j,k,3)*Grho(i,j,k)+Lvol(i,j,k,3)*Lrho(i,j  ,k)
           vol_l=Gvol(i,j-1,k,4)+Lvol(i,j-1,k,4); rho_l=Gvol(i,j-1,k,4)*Grho(i,j-1,k)+Lvol(i,j-1,k,4)*Lrho(i,j-1,k)
           RHOy (i,j,k)=(rho_l           +rho_r         )/(vol_l+vol_r)
	   betay(i,j,k)= 1.0d0/RHOy(i,j,k)
        end do
     end do
  end do

  do k=lo(3),hi(3)+1
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
	   ! Z face
           vol_r=Gvol(i,j,k,5)+Lvol(i,j,k,5); rho_r=Gvol(i,j,k,5)*Grho(i,j,k)+Lvol(i,j,k,5)*Lrho(i,j,k)
           vol_l=Gvol(i,j,k-1,6)+Lvol(i,j,k-1,6); rho_l=Gvol(i,j,k-1,6)*Grho(i,j,k-1)+Lvol(i,j,k-1,6)*Lrho(i,j,k-1)
           RHOz (i,j,k)=(rho_l           +rho_r         )/(vol_l+vol_r)
	   betaz(i,j,k)= 1.0d0/RHOz(i,j,k)
        end do
     end do
  end do

end subroutine interp_rho
 
end module mod_compute_interp_rho
